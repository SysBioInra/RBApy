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
        # WARNING: not storing document in a variable will result
        # in segmentation fault!
        document = self._load_document(input_file)
        model = document.model
        self._initialize_species(model, external_ids)
        self._initialize_reactions_and_enzymes(model)
        self._initialize_external_metabolites(model)
        self._initialize_transporter_information(model, cytosol_id)

    def _load_document(self, input_file):
        document = libsbml.readSBML(input_file)
        if document.getNumErrors() > 0:
            document.printErrors()
            raise UserWarning('Invalid SBML.')
        return document

    def _initialize_species(self, model, external_ids):
        if external_ids is None:
            external_ids = []
        self.species = rba.xml.ListOfSpecies()
        for spec in model.species:
            boundary = spec.getBoundaryCondition()
            if spec.getCompartment() in external_ids:
                boundary = True
            self.species.append(rba.xml.Species(spec.getId(), boundary))

    def _initialize_reactions_and_enzymes(self, model):
        reaction_list = extract_reactions(model)
        enzyme_list = self._extract_enzymes(model)
        self.reactions, self.enzymes \
            = duplicate_reactions_and_enzymes(reaction_list, enzyme_list)

    def _extract_enzymes(self, model):
        enzymes = FbcAnnotationParser(model).parse_enzymes()
        if not enzymes:
            enzymes = CobraNoteParser(model).read_notes()
        if not enzymes:
            print('Your SBML file does not contain fbc gene products nor uses '
                  'notes to define enzyme composition. Please comply with SBML'
                  ' requirements defined in the README and rerun script.')
            raise UserWarning('Invalid SBML.')
        return enzymes

    def _initialize_external_metabolites(self, model):
        self._tag_external_metabolites(model)
        self.external_metabolites = [m.id for m in self.species
                                     if m.boundary_condition]

    def _tag_external_metabolites(self, model):
        """Find metabolites not already tagged as external."""
        external_compartments = self._identify_external_compartments(model)
        for metabolite in self.species:
            compartment = model.getSpecies(metabolite.id).compartment
            if compartment in external_compartments:
                metabolite.boundary_condition = True

    def _identify_external_compartments(self, model):
        """
        External metabolites must meet the following conditions:
         (i) they participate in a sink/production reaction.
         (ii) all metabolites of their compartment meet condition (i).
        """
        sink_species = self._sink_species(model.reactions)
        result = set(c.id for c in model.compartments)
        for metabolite in model.species:
            if metabolite.id not in sink_species:
                result.discard(metabolite.compartment)
        return result

    def _sink_species(self, reactions):
        result = []
        for reaction in reactions:
            if (len(reaction.reactants) == 1 and
                    len(reaction.products) == 0):
                result.append(reaction.reactants[0].species)
            elif (len(reaction.products) == 1 and
                    len(reaction.reactants) == 0):
                result.append(reaction.products[0].species)
        return set(result)

    def _initialize_transporter_information(self, model, cytosol_id):
        self.imported_metabolites = self._transport_reactions(
            model, cytosol_id
        )
        self.has_membrane_enzyme = find_membrane_reactions(self.reactions)

    def _transport_reactions(self, model, cytosol_id):
        """
        Identify all transport reactions in the SBML file.

        They meet the following conditions:
        - one of the reactants has the same prefix (e.g. M_glc) as one of the
        external metabolites. This product should not be in the cytosol.
        - one of the products is in the cytosol.
        """
        external_prefixes = set(
            self._prefix(m) for m in self.external_metabolites
        )
        result = {}
        for reaction in self.reactions:
            if self._has_cytosolic_product(reaction, cytosol_id):
                transported = self._noncytosolic_external_reactants(
                    reaction, cytosol_id, external_prefixes
                )
                if transported:
                    result[reaction.id] = transported
        return result

    def _prefix(self, metabolite_id):
        return metabolite_id.rsplit('_', 1)[0]

    def _has_cytosolic_product(self, reaction, cytosol_id):
        return any(self._suffix(p.species) == cytosol_id
                   for p in reaction.products)

    def _suffix(self, metabolite_id):
        return metabolite_id.rsplit('_', 1)[1]

    def _noncytosolic_external_reactants(self, reaction, cytosol_id,
                                         external_prefixes):
        result = []
        for reactant in reaction.reactants:
            prefix, cpt = reactant.species.rsplit('_', 1)
            if cpt != cytosol_id and prefix in external_prefixes:
                result.append(reactant.species)
        return result


def extract_reactions(model):
    """
    Parse reactions.

    Returns
    -------
    rba.xml.ListOfReactions
        Reactions in RBA format.

    """
    reactions = rba.xml.ListOfReactions()
    for reaction in model.reactions:
        new_reaction = rba.xml.Reaction(reaction.id, reaction.reversible)
        for reactant in reaction.reactants:
            new_reaction.reactants.append(
                rba.xml.SpeciesReference(reactant.species,
                                         reactant.stoichiometry)
            )
        for product in reaction.products:
            new_reaction.products.append(
                rba.xml.SpeciesReference(product.species,
                                         product.stoichiometry)
            )
        reactions.append(new_reaction)
    return reactions


def duplicate_reactions_and_enzymes(reaction_list, enzyme_list):
    """
    Duplicate reactions catalyzed by several enzymes.

    Parameters
    ----------
    reactions: rba.xml.ListOfReactions
        List of reactions.

    enzymes: list of lists of lists
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


class FbcAnnotationParser(object):
    """Parse fbc annotation to gather enzyme compositions."""
    def __init__(self, model):
        self._model = model
        self._fbc = model.getPlugin('fbc')

    def parse_enzymes(self):
        if not self._fbc:
            return []
        self._initialize_gene_id_to_name_map()
        return [self._enzyme_composition(r) for r in self._model.reactions]

    def _initialize_gene_id_to_name_map(self):
        self._gene_names = {}
        for gene_product in self._fbc.getListOfGeneProducts():
            self._gene_names[gene_product.id] = gene_product.label

    def _enzyme_composition(self, reaction):
        gp_association = reaction.getPlugin('fbc') \
                                 .getGeneProductAssociation()
        if gp_association:
            return self._read_fbc_association(gp_association)
        else:
            return [[]]

    def _read_fbc_association(self, gp_association):
        """We assume that relations are always 'or's of 'and's."""
        association = gp_association.getAssociation()
        if association.isGeneProductRef():
            gene_id = association.getGeneProduct()
            return [[self._gene_names[gene_id]]]
        elif association.isFbcOr():
            return [self._read_fbc_association_components(a)
                    for a in association.getListOfAssociations()]
        elif association.isFbcAnd():
            result = []
            for assoc in association.getListOfAssociations():
                result += self._read_fbc_association_components(assoc)
            return [result]
        else:
            print('Invalid association.')
            raise UserWarning('Invalid SBML.')

    def _read_fbc_association_components(self, association):
        """We assume that relations are always 'and's."""
        if association.isGeneProductRef():
            gene_id = association.getGeneProduct()
            return [self._gene_names[gene_id]]
        elif association.isFbcAnd():
            result = []
            for assoc in association.getListOfAssociations():
                result += self._read_fbc_association_components(assoc)
            return result
        else:
            print('Invalid association (we only support ors of ands)')
            raise UserWarning('Invalid SBML.')


class CobraNoteParser(object):
    def __init__(self, model):
        self._model = model

    def read_notes(self):
        reactions = self._model.reactions
        result = []
        for reaction in self._model.reactions:
            if not reaction.notes:
                return []
            result.append(self._parse_note(reaction.notes))
        return result

    def _parse_note(self, note):
        result = []
        for ga in self._gene_associations(note):
            composition = self._parse_gene_association(ga)
            if composition:
                result.append(composition)
        return result

    def _gene_associations(self, note):
        # fields may be encapsulated in a <html> tag (or equivalent)
        note = self._remove_html_tag(note)
        return (note.getChild(i).getChild(0).toString()
                for i in range(notes.getNumChildren()))

    def _remove_html_tag(self, note):
        if (note.getNumChildren() == 1
                and note.getChild(0).getName() != "p"):
            return note.getChild(0)
        return note

    def _parse_gene_association(self, text):
        """We assume that relations are always 'or's of 'and's."""
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
