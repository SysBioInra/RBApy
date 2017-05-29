
import libsbml
import copy
import itertools

from ..rba_xml import *

class SBMLFilter:
    """
    Class used to filter RBA-relevant SBML data.
    """
    def __init__(self, input_file, cytosol_id = 'c', external_ids = []):
        """
        Constructor from file path.

        :param input: Path to input file.
        :type input: String.
        """
        # load SBML file
        input_document = libsbml.readSBML(input_file)
        if input_document.getNumErrors() > 0:
            input_document.printErrors()
            raise UserWarning('Invalid SBML.')

        # store metabolite info
        self.species = ListOfSpecies()
        for s in input_document.getModel().species:
            boundary_condition = s.getBoundaryCondition()
            if s.getCompartment() in external_ids:
                boundary_condition = True
            self.species.append(Species(s.getId(), boundary_condition))

        # extract enzymes associated to reactions
        self.enzymes = []
        self.reactions = ListOfReactions()
        self._extract_enzymes_and_reactions(input_document)

        # find boundary conditions
        self._find_boundary_condition(input_document)
        self.external_metabolites = [m.id for m in self.species \
                                     if m.boundary_condition]

        # identify membrane and transport reactions
        self.imported_metabolites = []
        self._find_transport_reactions(cytosol_id)
        self._find_membrane_reactions()
                    
    def _extract_enzymes_and_reactions(self, input_document):
        """
        Parse annotation containing enzyme components.
        """
        self.enzymes = []
        self.reactions = ListOfReactions()
        # try to read fbc notes or read old-fashioned notes
        enzyme_list = self._read_fbc_annotation(input_document)
        if not(enzyme_list):
            enzyme_list =  self._read_notes(input_document)
        if enzyme_list:
            reactions = input_document.getModel().reactions
            for index in range(len(reactions)):
                # create reaction in RBA objects
                r = reactions[index]
                new_reaction = Reaction(r.getId(), r.getReversible());
                for reactant in r.reactants:
                    new_reaction.reactants.append \
                        (SpeciesReference(reactant.getSpecies(),
                                          reactant.getStoichiometry()))
                for product in r.products:
                    new_reaction.products.append \
                        (SpeciesReference(product.getSpecies(),
                                          product.getStoichiometry()))
                # duplicate reactions having multiple enzymes
                suffix = 0
                for enzyme in enzyme_list[index]:
                    suffix += 1
                    reaction = copy.copy(new_reaction)
                    reaction.id += '_' + str(suffix)
                    self.enzymes.append(enzyme)
                    self.reactions.append(reaction)
        else:
            print('Your SBML file does not contain fbc gene products nor uses '
                  'notes to define enzyme composition. Please comply with SBML '
                  'requirements defined in the README and rerun script.')
            raise UserWarning('Invalid SBML.')

    def _find_boundary_condition(self, input_document):
        """
        Go through reactions and identify all external metabolites that were
        not tagged as having boundaryCondition=True. 
        External metabolites must meet the following conditions:
         (i) they participate in a sink/production reaction.
         (ii) all metabolites of their compartment meet condition (i).
        """
        # identify species participating in a sink reaction
        sink_species = []
        for r in input_document.getModel().reactions:
            if len(r.reactants) + len(r.products) > 1: continue
            if len(r.reactants) == 1:
                sink_species.append(r.reactants[0].getSpecies())
            else:
                sink_species.append(r.products[0].getSpecies())
        # find external compartments
        external_compartments = [c.getId() for c in input_document.getModel().compartments]
        for metabolite in input_document.getModel().species:
            if metabolite.getId() not in sink_species:
                try:
                    external_compartments.remove(metabolite.getCompartment())
                except ValueError:
                    pass
        # tag all metabolites belonging to external compartments as boundary
        for metabolite in self.species:
            compartment = input_document.getModel().getSpecies(metabolite.id).compartment
            if compartment in external_compartments:
                metabolite.boundary_condition = True

    def _find_membrane_reactions(self):
        """
        Identify all reactions whose enzyme is at least partly in the membrane.
        """
        self.has_membrane_enzyme = {}
        for r in self.reactions:
            compartments = [m.species.rsplit('_',1)[1]
                            for m in itertools.chain(r.reactants, r.products)]
            self.has_membrane_enzyme[r.id] \
                = any(c != compartments[0] for c in compartments[1:])

    def _find_transport_reactions(self, cytosol_id):
        """
        Identify all transport reactions in the SBML file. They meet the
        following conditions:
        - one of the products has the same prefix (e.g. M_glc) as one of the
        external metabolites. This product should not be in the cytosol.
        - one of the reactants is in the cytosol.
        """
        external_prefixes \
            = [m.rsplit('_',1)[0] for m in self.external_metabolites]
        self.imported_metabolites = {}
        for r in self.reactions:
            transported = []
            # check that one of the products is in the cytosol
            comps = [p.species.rsplit('_',1)[1] for p in r.products]
            if all(comp != cytosol_id for comp in comps): continue
            # look if one of the reactant has the prefix of an external
            # metabolite and is NOT in the cytosol
            transported = []
            for m in r.reactants:
                [prefix, comp] = m.species.rsplit('_',1)
                if comp != cytosol_id and prefix in external_prefixes:
                    transported.append(m.species)
            if transported:
                self.imported_metabolites[r.id] = transported
            
    def _read_notes(self, input_document):
        """
        Parse old-fashioned notes containing enzyme components.
        """
        reactions = input_document.getModel().reactions
        enzymes = []
        for r in reactions:
            notes = r.getNotes()
            # check that a note is indeed available
            if not(notes): return None
            # fields may be encapsulated in a <html> tag (or equivalent)
            if notes.getNumChildren() == 1 \
               and notes.getChild(0).getName() != "p":
                notes = notes.getChild(0)
            for i in range(notes.getNumChildren()):
                text = notes.getChild(i).getChild(0).toString()
                enzyme_composition = read_gene_association(text)
                if enzyme_composition:
                    enzymes.append(enzyme_composition)
        return enzymes

    def _read_fbc_annotation(self, input_document):
        """
        Parse fbc annotation to gather enzyme compositions.
        """
        # get fbc annotation (if available)
        fbc = input_document.getModel().getPlugin('fbc')
        if not(fbc): return None

        # get gene id - gene name association
        gene_names = {}
        for gene_product in fbc.gene_products:
            gene_names[gene_product.getId()] = gene_product.getLabel()

        # gather enzyme composition
        enzyme_list = []
        reactions = input_document.getModel().reactions
        for r in reactions:
            # get fbc:geneProductAssociation
            gp_association = r.getPlugin('fbc').getGeneProductAssociation()
            if gp_association:
                enzymes = read_fbc_association(gp_association, gene_names)
            else:
                # no gene association: we assume this reaction is
                # spontaneous
                enzymes = [[]]
            enzyme_list.append(enzymes)
        
        return enzyme_list

def read_gene_association(text):
    """
    Parse GENE_ASSOCIATION field from string representing a note field.
    (for this version, we assume that relations are always 'or's of 'and's)

    :params text: Note field containing GENE_ASSOCIATION.
    :type text: string
    """
    tags = text.split(':',1)
    if len(tags) != 2:
        print('Invalid note field: ' + text)
        return None
    if tags[0] != "GENE_ASSOCIATION":
        return None
    else:
        enzyme_set = tags[1]
        if len(enzyme_set) == 0: return []
        
        # field is not standard: we try to standardize a little.
        # remove parentheses
        enzyme_set = ''.join(c for c in enzyme_set if c not in '()')
        
        enzymes = enzyme_set.split(' or ')
        compositions = []
        for enzyme in enzymes:
            compositions.append(list(map(str.strip, enzyme.split(' and '))))
        return compositions

def read_fbc_association(gp_association, gene_names = None):
    """
    Parse fbc:geneProductAssociation and return the list of
    proteins composing enzyme.
    (for this version, we assume that relations are always 'or's of 'and's)

    :parameter gp_association: gene product assocation object to parse.
    :parameter gene_names: dictionary used to replace gene ids by their name.
    :type gp_association: FbcGeneProductAssociation instance
    :type gene_names: dictionary
    """
    association = gp_association.getAssociation()
    if association.isGeneProductRef():
        gene_id = association.getGeneProduct()
        if gene_names:
            return [[gene_names[gene_id]]]
        else:
            return [[gene_id]]
    elif association.isFbcOr():
        return [read_fbc_association_components(a, gene_names) \
                for a in association.getListOfAssociations()]
    elif association.isFbcAnd():
        result = []
        for a in association.getListOfAssociations():
            result += read_fbc_association_components(a, gene_names)
        return [result]
    else:
        print('Invalid association.')
        raise UserWarning('Invalid SBML.')

def read_fbc_association_components(association, gene_names = None):
    """
    Parse fbc:Association and return the list of gene ids it contains.
    (for this version, we assume that relations are always 'and's)

    :parameter association: assocation object to parse.
    :parameter gene_names: dictionary used to replace gene ids by their name.
    :type association: FbcAssociation instance
    :type gene_names: dictionary
    """
    if association.isGeneProductRef():
        gene_id = association.getGeneProduct()
        if gene_names:
            return [gene_names[gene_id]]
        else:
            return [gene_id]
    elif association.isFbcAnd():
        result = []
        for a in association.getListOfAssociations():
            result += read_fbc_association_components(a, gene_names)
        return result
    else:
        print('Invalid association (well not really but I was hoping it would '
              'be ors of ands :/')
        raise UserWarning('Invalid SBML.')
        
