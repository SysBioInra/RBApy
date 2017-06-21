"""
Module defining SBMLFilter class.
"""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import copy
import itertools
import libsbml

# local imports
from rba import rba_xml

class SBMLFilter(object):
    """
    Class used to filter RBA-relevant SBML data.

    Attributes:
        species: rba_xml.ListOfSpecies containing SBML species.
        enzymes: list where each element represents an enzyme as a list of
            protein identifiers.
        reactions: rba_xml.ListOfReaction containing SBML reactions.
        external_metabolites: list of SBML identifiers of external metabolites.
        imported_metabolites: list of SBML identifiers of imported metabolites.
        has_membrane_enzyme: dict mapping reaction id to boolean value
            indicating whether reaction occurs in the membrane.
    """

    def __init__(self, input_file, cytosol_id='c', external_ids=None):
        """
        Constructor from file.

        Args:
            input: Path to input file.
            cytosol_id: identifier of cytosol in the SBML file.
            external_ids: identifiers of external compartments in SBML file.
        """
        # load SBML file
        input_document = libsbml.readSBML(input_file)
        if input_document.getNumErrors() > 0:
            input_document.printErrors()
            raise UserWarning('Invalid SBML.')

        # store metabolite info
        if external_ids is None:
            external_ids = []
        self.species = rba_xml.ListOfSpecies()
        for spec in input_document.getModel().species:
            boundary = spec.getBoundaryCondition()
            if spec.getCompartment() in external_ids:
                boundary = True
            self.species.append(rba_xml.Species(spec.getId(), boundary))

        # extract enzymes associated to reactions
        self.enzymes, self.reactions \
            = self._extract_enzymes_and_reactions(input_document)

        # find external metabolites
        self._set_boundary_condition(input_document)
        self.external_metabolites = [m.id for m in self.species
                                     if m.boundary_condition]

        # identify membrane and transport reactions
        self.imported_metabolites = self._find_transport_reactions(cytosol_id)
        self.has_membrane_enzyme = self._find_membrane_reactions()

    def _extract_enzymes_and_reactions(self, input_document):
        """
        Parse annotation containing enzyme components.
        """
        enzymes = []
        reactions = rba_xml.ListOfReactions()
        # try to read fbc notes, else read old-fashioned notes
        enzyme_list = self._read_fbc_annotation(input_document)
        if not enzyme_list:
            enzyme_list = self._read_notes(input_document)
        if enzyme_list:
            sbml_reactions = input_document.getModel().reactions
            for index, reaction in enumerate(sbml_reactions):
                # create reaction in RBA objects
                new_reaction = rba_xml.Reaction(reaction.getId(),
                                                reaction.getReversible())
                for reactant in reaction.reactants:
                    sr = rba_xml.SpeciesReference(reactant.getSpecies(),
                                                  reactant.getStoichiometry())
                    new_reaction.reactants.append(sr)
                for product in reaction.products:
                    sr = rba_xml.SpeciesReference(product.getSpecies(),
                                                  product.getStoichiometry())
                    new_reaction.products.append(sr)
                # duplicate reactions having multiple enzymes
                suffix = 0
                for enzyme in enzyme_list[index]:
                    suffix += 1
                    clone = copy.copy(new_reaction)
                    clone.id += '_' + str(suffix)
                    enzymes.append(enzyme)
                    reactions.append(clone)
        else:
            print('Your SBML file does not contain fbc gene products nor uses '
                  'notes to define enzyme composition. Please comply with SBML '
                  'requirements defined in the README and rerun script.')
            raise UserWarning('Invalid SBML.')
        return enzymes, reactions

    def _set_boundary_condition(self, input_document):
        """
        Go through reactions and identify external metabolites.

        We try to find metabolites whose boundaryCondition is not already set.
        External metabolites must meet the following conditions:
         (i) they participate in a sink/production reaction.
         (ii) all metabolites of their compartment meet condition (i).
        """
        # identify species participating in a sink reaction
        sink_species = []
        for reaction in input_document.getModel().reactions:
            if len(reaction.reactants) + len(reaction.products) > 1:
                continue
            if len(reaction.reactants) == 1:
                sink_species.append(reaction.reactants[0].getSpecies())
            else:
                sink_species.append(reaction.products[0].getSpecies())
        # find external compartments
        external_compartments \
            = [c.getId() for c in input_document.getModel().compartments]
        for metabolite in input_document.getModel().species:
            if metabolite.getId() not in sink_species:
                try:
                    external_compartments.remove(metabolite.getCompartment())
                except ValueError:
                    pass
        # tag all metabolites belonging to external compartments as boundary
        for metabolite in self.species:
            compartment = (input_document.getModel()
                           .getSpecies(metabolite.id).compartment)
            if compartment in external_compartments:
                metabolite.boundary_condition = True

    def _find_membrane_reactions(self):
        """
        Identify all reactions whose enzyme is at least partly in the membrane.
        """
        result = {}
        for reaction in self.reactions:
            compartments = [m.species.rsplit('_', 1)[1]
                            for m in itertools.chain(reaction.reactants,
                                                     reaction.products)]
            result[reaction.id] \
                = any(c != compartments[0] for c in compartments[1:])
        return result

    def _find_transport_reactions(self, cytosol_id):
        """
        Identify all transport reactions in the SBML file.

        They meet the following conditions:
        - one of the products has the same prefix (e.g. M_glc) as one of the
        external metabolites. This product should not be in the cytosol.
        - one of the reactants is in the cytosol.
        """
        external_prefixes \
            = [m.rsplit('_', 1)[0] for m in self.external_metabolites]
        imported_metabolites = {}
        for reaction in self.reactions:
            transported = []
            # check that one of the products is in the cytosol
            comps = [p.species.rsplit('_', 1)[1] for p in reaction.products]
            if all(comp != cytosol_id for comp in comps):
                continue
            # look if one of the reactant has the prefix of an external
            # metabolite and is NOT in the cytosol
            transported = []
            for reactant in reaction.reactants:
                [prefix, comp] = reactant.species.rsplit('_', 1)
                if comp != cytosol_id and prefix in external_prefixes:
                    transported.append(reactant.species)
            if transported:
                imported_metabolites[reaction.id] = transported
        return imported_metabolites

    @staticmethod
    def _read_notes(input_document):
        """
        Parse old-fashioned notes containing enzyme components.
        """
        reactions = input_document.getModel().reactions
        enzymes = []
        for reaction in reactions:
            notes = reaction.getNotes()
            # check that a note is indeed available
            if not notes: return None
            # fields may be encapsulated in a <html> tag (or equivalent)
            if (notes.getNumChildren() == 1
                and notes.getChild(0).getName() != "p"):
                notes = notes.getChild(0)
            for i in range(notes.getNumChildren()):
                text = notes.getChild(i).getChild(0).toString()
                enzyme_composition = read_gene_association(text)
                if enzyme_composition:
                    enzymes.append(enzyme_composition)
        return enzymes

    @staticmethod
    def _read_fbc_annotation(input_document):
        """
        Parse fbc annotation to gather enzyme compositions.
        """
        # get fbc annotation (if available)
        fbc = input_document.getModel().getPlugin('fbc')
        if not fbc: return None

        # get gene id - gene name association
        gene_names = {}
        for gene_product in fbc.getListOfGeneProducts():
            gene_names[gene_product.getId()] = gene_product.getLabel()

        # gather enzyme composition
        enzyme_list = []
        reactions = input_document.getModel().reactions
        for reaction in reactions:
            # get fbc:geneProductAssociation
            gp_association = reaction.getPlugin('fbc') \
                                     .getGeneProductAssociation()
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

    For this version, we assume that relations are always 'or's of 'and's.

    Args:
        text: Note field containing GENE_ASSOCIATION.

    Returns:
        List of enzymes extracted. Every enzyme is represented as a list of
        protein identifiers.
    """
    tags = text.split(':', 1)
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
        # split enzymes
        enzymes = enzyme_set.split(' or ')
        compositions = []
        for enzyme in enzymes:
            compositions.append([e.strip() for e in enzyme.split(' and ')])
        return compositions

def read_fbc_association(gp_association, gene_names=None):
    """
    Parse fbc:geneProductAssociation.

    For this version, we assume that relations are always 'or's of 'and's.

    Args:
        gp_association: standard fbc gene product assocation object.
        gene_names: dictionary used to replace gene ids by their name.

    Returns:
        List of enzymes extracted. Every enzyme is represented as a list of
        gene identifiers.
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
        for assoc in association.getListOfAssociations():
            result += read_fbc_association_components(assoc, gene_names)
        return [result]
    else:
        print('Invalid association.')
        raise UserWarning('Invalid SBML.')

def read_fbc_association_components(association, gene_names=None):
    """
    Parse fbc:Association and return the list of gene names it contains.

    For this version, we assume that relations are always 'and's.

    Args:
        association: fbc assocation object.
        gene_names: dictionary used to replace gene ids by their name.

    Returns:
        List of gene names.
    """
    if association.isGeneProductRef():
        gene_id = association.getGeneProduct()
        if gene_names:
            return [gene_names[gene_id]]
        else:
            return [gene_id]
    elif association.isFbcAnd():
        result = []
        for assoc in association.getListOfAssociations():
            result += read_fbc_association_components(assoc, gene_names)
        return result
    else:
        print('Invalid association (well not really but I was hoping it would '
              'be ors of ands :/')
        raise UserWarning('Invalid SBML.')
