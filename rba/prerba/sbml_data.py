"""Module defining SbmlData class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import copy
import itertools
import libsbml

# local imports
import rba.xml


class SbmlData(object):
    """
    Class used to parse RBA-relevant SBML data.

    Attributes
    ----------
    species: rba.xml.ListOfSpecies
        SBML species.
    enzymes: list
        Enzymes as a list of protein identifiers.
    reactions: rba.xml.ListOfReaction
        SBML reactions.
    external_metabolites: list
        SBML identifiers of external metabolites.
    imported_metabolites: list
        SBML identifiers of imported metabolites.
    has_membrane_enzyme: dict
        Keys are reaction ids and values booleans
        indicating whether reaction occurs in the membrane.

    """

    def __init__(self, input_file, cytosol_id='c', external_ids=None):
        """
        Build from file.

        Parameters
        ----------
        input: str
            Path to input file.
        cytosol_id: str
            identifier of cytosol in the SBML file.
        external_ids: list of str
            identifiers of external compartments in SBML file.

        """
        document = self._load_document(input_file)
        self._initialize_species(document, external_ids)
        self._initialize_reactions_and_enzymes(document)
        self._initialize_external_metabolites(document)
        self._initialize_transporter_information(document, cytosol_id)

    def _load_document(self, input_file):
        document = libsbml.readSBML(input_file)
        if document.getNumErrors() > 0:
            document.printErrors()
            raise UserWarning('Invalid SBML.')
        return document

    def _initialize_species(self, input_document, external_ids):
        if external_ids is None:
            external_ids = []
        self.species = rba.xml.ListOfSpecies()
        for spec in input_document.getModel().species:
            boundary = spec.getBoundaryCondition()
            if spec.getCompartment() in external_ids:
                boundary = True
            self.species.append(rba.xml.Species(spec.getId(), boundary))

    def _initialize_reactions_and_enzymes(self, input_document):
        reaction_list = extract_reactions(input_document)
        enzyme_list = extract_enzymes(input_document)
        self.reactions, self.enzymes \
            = duplicate_reactions_and_enzymes(reaction_list, enzyme_list)

    def _initialize_external_metabolites(self, input_document):
        self._set_boundary_condition(input_document)
        self.external_metabolites = [m.id for m in self.species
                                     if m.boundary_condition]

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

    def _initialize_transporter_information(self, input_document, cytosol_id):
        self.imported_metabolites = self._find_transport_reactions(cytosol_id)
        self.has_membrane_enzyme = find_membrane_reactions(self.reactions)

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


def extract_reactions(input_document):
    """
    Parse reactions.

    Parameters
    ----------
    input_document: SBML document
        Valid SBML document.

    Returns
    -------
    rba.xml.ListOfReactions
        Reactions in RBA format.

    """
    reactions = rba.xml.ListOfReactions()
    for reaction in input_document.model.reactions:
        # create reaction in RBA objects
        new_reaction = rba.xml.Reaction(reaction.getId(),
                                        reaction.getReversible())
        for reactant in reaction.reactants:
            sr = rba.xml.SpeciesReference(reactant.getSpecies(),
                                          reactant.getStoichiometry())
            new_reaction.reactants.append(sr)
        for product in reaction.products:
            sr = rba.xml.SpeciesReference(product.getSpecies(),
                                          product.getStoichiometry())
            new_reaction.products.append(sr)
        reactions.append(new_reaction)
    return reactions


def duplicate_reactions_and_enzymes(reaction_list, enzyme_list):
    """
    Duplicate reactions catalyzed by several enzymes.

    Parameters
    ----------
    reaction_list: rba.xml.ListOfReactions
        List of reactions.

    enzyme_list: list of lists of lists
        List of same length as reactions. Every reaction should be represented
        as the list of enzyme that catalyzes it. Every enzyme should be
        reprensented as the list of proteins that compose it.

    Returns
    -------
    (reactions, enzymes) couple with
    reactions: rba.xml.ListOfReactions
        Updated list of reactions where reactions with multiple enzymes have
        been duplicated.
    enzymes: list of lists
        Updated list of enzymes. Now every reaction is represented by exactly
        one enzyme.

    """
    enzymes = []
    reactions = rba.xml.ListOfReactions()
    for reaction, reaction_enzymes in zip(reaction_list, enzyme_list):
        # duplicate reactions having multiple enzymes
        suffix = 0
        for enzyme in reaction_enzymes:
            suffix += 1
            clone = copy.copy(reaction)
            if suffix > 1:
                clone.id += '_' + str(suffix)
            enzymes.append(enzyme)
            reactions.append(clone)
    return reactions, enzymes


def find_reaction_location(reaction):
    """
    Find reactions with reactants and products in a single compartment.

    Parameters
    ----------
    reactions: rba.xml.Reaction
        Reaction.

    Returns
    -------
    str
        Compartment containing reaction or empty string if reactants and
        products were in more than one compartment.

    """
    compartments = [m.species.rsplit('_', 1)[1]
                    for m in itertools.chain(reaction.reactants,
                                             reaction.products)]
    return all(c == compartments[0] for c in compartments[1:])


def enzymatic_protein_location(enzyme_comp, reactions):
    """
    Guess location of enzymatic proteins.

    For every protein we record the location of the reactions it catalyzes.
    If locations match, we assume that the protein is contained
    in this compartment.

    Parameters
    ----------
    enzyme_comp: list of lists
        List of enzymes where an enzyme is a list of protein ids.
    reactions: rba.xml.ListOfReactions
        List of reactions.

    Returns
    -------
    dict
        Dictionary containing locations predicted. Key are protein ids and
        values are compartment names.

    """
    protein_location_list = {}
    for enzyme_composition, reaction in zip(enzyme_comp, reactions):
        location = sbml_filter.find_reaction_location(reaction)
        if location:
            for prot in enzyme_comp:
                protein_location_list.setdefault(prot, []).append(location)
    result = {}
    for protein, location_list in protein_location_list:
        if (len(location_list) == 1 or
                all(loc == location_list[0] for loc in location_list[1:])):
            result[protein] = location_list[0]
    return result


def find_membrane_reactions(reactions):
    """
    Identify reactions with membrane enzymes.

    Parameters
    ----------
    reactions: rba.xml.ListOfReactions
        List of reactions.

    Returns
    -------
    dict
        Dictionary where keys are reactions. Value associated to key
        is true if enzyme associated to reaction is a membrane enzyme,
        false otherwise.

    """
    result = {}
    for reaction in reactions:
        compartments = [m.species.rsplit('_', 1)[1]
                        for m in itertools.chain(reaction.reactants,
                                                 reaction.products)]
        result[reaction.id] \
            = any(c != compartments[0] for c in compartments[1:])
    return result


def extract_enzymes(input_document):
    """
    Parse annotation containing enzyme components.

    Parameters
    ----------
    input_document: SBML document
        Valid SBML document.

    Returns
    -------
    list of lists of lists
        For every reaction in the model, return a list of enzymes that
        catalyze this reaction. An enzyme is represented by a list of
        proteins that composes it.

    """
    # try to read fbc notes, else read old-fashioned notes
    enzymes = read_fbc_annotation(input_document)
    if not enzymes:
        enzymes = read_notes(input_document)
    if not enzymes:
        print('Your SBML file does not contain fbc gene products nor uses '
              'notes to define enzyme composition. Please comply with SBML'
              ' requirements defined in the README and rerun script.')
        raise UserWarning('Invalid SBML.')
    return enzymes


def read_notes(input_document):
    """
    Parse old-fashioned notes containing enzyme components.

    Parameters
    ----------
    input_document: SBML document
        Valid SBML document.

    Returns
    -------
    list of lists of lists
        For every reaction in the model, return a list of enzymes that
        catalyze this reaction. An enzyme is represented by a list of
        proteins that composes it.

    """
    reactions = input_document.getModel().reactions
    enzymes = []
    for reaction in reactions:
        notes = reaction.getNotes()
        # check that a note is indeed available
        if not notes:
            return None
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


def read_fbc_annotation(input_document):
    """
    Parse fbc annotation to gather enzyme compositions.

    Parameters
    ----------
    input_document: SBML document
        Valid SBML document.

    Returns
    -------
    list of lists of lists
        For every reaction in the model, return a list of enzymes that
        catalyze this reaction. An enzyme is represented by a list of
        proteins that composes it.

    """
    # get fbc annotation (if available)
    fbc = input_document.getModel().getPlugin('fbc')
    if not fbc:
        return None

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

    Parameters
    ----------
    text: string
        Note field containing GENE_ASSOCIATION.

    Returns
    -------
    list of lists
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
        if len(enzyme_set) == 0:
            return []
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

    Parameters
    ----------
    gp_association: fbc:geneProductAssociation object
        association to parse.
    gene_names: dict
        dictionary used to replace gene ids by their name.

    Returns
    -------
    list of lists
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
        return [read_fbc_association_components(a, gene_names)
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

    Parameters
    ----------
    association: fbc:Assocation object
        association to parse.
    gene_names: object
        dictionary used to replace gene ids by their name.

    Returns
    -------
    list
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
        print('Invalid association (we only support ors of ands)')
        raise UserWarning('Invalid SBML.')
