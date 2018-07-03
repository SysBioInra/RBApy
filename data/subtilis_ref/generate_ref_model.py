
from lxml import etree
import libsbml
from itertools import chain
import sys
from os.path import join

sys.path.append(join(sys.path[0], '../..'))

import rba  # noqa
from rba import rba_xml  # noqa


def main():
    output = '.'
    ModelConverter('old_data/').model.write_files('.')


class ModelConverter(object):
    def __init__(self, input_dir):
        self._input_dir = input_dir
        self.model = rba.RbaModel()
        self._read_metabolism()
        self.model.proteins = self._read_macromolecules('proteins.xml',
                                                        'protein')
        self.model.rnas = self._read_macromolecules('rnas.xml', 'rna')
        self.model.dna = self._read_macromolecules('dna.xml', 'dna')
        self._read_parameters()
        self._read_enzymes()
        self._read_processes()
        self.model.medium = self.model.read_medium('medium.tsv')
        self.adapt_metabolite_names()
        self.add_matp_process()

    def _read_metabolism(self):
        sbml = libsbml.readSBML(self._input_dir + 'metabolism.xml')
        metabolism = self.model.metabolism
        for c in sbml.model.compartments:
            metabolism.compartments.append(rba_xml.Compartment(c.id))
        for s in sbml.model.species:
            metabolism.species.append(rba_xml.Species(s.id,
                                                      s.boundary_condition))
        for r in sbml.model.reactions:
            new_reaction = rba_xml.Reaction(r.id, r.reversible)
            metabolism.reactions.append(new_reaction)
            for sr in r.reactants:
                new_reaction.reactants.append(
                    rba_xml.SpeciesReference(sr.species, sr.stoichiometry)
                    )
            for sr in r.products:
                new_reaction.products.append(
                    rba_xml.SpeciesReference(sr.species, sr.stoichiometry)
                    )

    def _read_macromolecules(self, filename, tag):
        root = self._xml_root(filename)
        result = rba_xml.RbaMacromolecules()
        result.components = rba_xml.ListOfComponents.from_xml_node(
            rba_xml.get_unique_child(root, 'listOfComponents')
        )
        for m in root.find('listOfSpecies').iterfind(tag):
            result.macromolecules.append(
                rba_xml.Macromolecule.from_xml_node(m)
                )
        return result

    def _xml_root(self, filename):
        return etree.ElementTree(file=self._input_dir + filename).getroot()

    def _read_parameters(self):
        self._read_functions()
        self._read_densities()

    def _read_functions(self):
        root = self._xml_root('parameters.xml')
        func = rba_xml.get_unique_child(root, 'listOfFunctions')
        self.model.parameters.functions = rba_xml.ListOfFunctions.from_xml_node(func)

    def _read_densities(self):
        root = self._xml_root('parameters.xml')
        for d in root.iterfind('listOfMaximalDensities/maximalDensity'):
            target = rba_xml.TargetDensity(d.get('compartment'))
            target.upper_bound = self._parse_value(
                d, target.compartment + '_density'
                )
            self.model.density.target_densities.append(target)

    def _parse_value(self, node, id_):
        """
        Parse old-fashioned value node.

        If single function reference, return it, else create function
        (or aggregate) with proposed id and return id.
        """
        value = node.get('value')
        if value is not None:
            try:
                new_fn = rba_xml.Function(id_, 'constant',
                                          {'CONSTANT': float(value)})
            except ValueError:
                return value
            self.model.parameters.functions.append(new_fn)
            return id_
        fn_refs = node.findall('functionReference')
        if len(fn_refs) == 1:
            return fn_refs[0].get('function')
        new_agg = rba_xml.Aggregate(id_, 'multiplication')
        for ref in fn_refs:
            new_ref = rba_xml.FunctionReference(ref.get('function'))
            new_agg.function_references.append(new_ref)
        self.model.parameters.aggregates.append(new_agg)
        return id_

    def _read_enzymes(self):
        EnzymeConverter(self._xml_root('enzymes.xml'), self.model)

    def _read_processes(self):
        self.read_process_machineries()
        self.read_targets()
        self.read_component_maps()

    def _read_process_machineries(self):
        root = self._xml_root('processes.xml')
        for p in root.iterfind('listOfProcesses/process'):
            new_p = rba_xml.Process(p.get('id'), p.get('name'))
            self.model.processes.processes.append(new_p)
            capacity = p.find('capacityConstraint')
            if capacity is not None:
                new_p.machinery.machinery_composition \
                    = rba_xml.MachineryComposition.from_xml_node(
                        capacity.find('machineryComposition')
                    )
                new_p.machinery.capacity.value = self.parse_value(
                    capacity.find('capacity'), p.get('id') + '_capacity'
                )
            for prod in p.iterfind('operatingCosts/production'):
                new_p.processings.productions.append(
                    self._create_processing(prod)
                )
            for deg in p.iterfind('operatingCosts/degradation'):
                new_p.processings.degradations.append(
                    rba_xml.Operation.from_xml_node(deg)
                    )

    def _create_processing(self, node):
        result = rba_xml.Processing.from_xml_node(node)
        if result.set == 'protein':
            set_ = self.model.proteins
        elif result.set == 'rna':
            set_ = self.model.rnas
        elif resut.set == 'dna':
            set_ = self.model.dna
        for id_ in (m.id for m in set_.macromolecules):
            result.inputs.append(rba_xml.SpeciesReference(id_, 1))
        return result

    def _read_component_maps(self):
        for map_ in root.iterfind('listOfComponentMaps/componentMap'):
            new_map = rba_xml.ComponentMap(map_.get('id'))
            self.model.processes.component_maps.append(new_m)
            node = map_.find('constantCost')
            if node is not None:
                new_m.constant_cost = rba_xml.ConstantCost.from_xml_node(node)
            for c in map_.iterfind('cost'):
                new_m.costs.append(rba_xml.Cost.from_xml_node(c))

    def _read_targets(self):
        root = self._xml_root('processes.xml')
        for process in root.iterfind('listOfProcesses/process'):
            new_group = rba_xml.TargetGroup(process.get('id') + '_targets')
            for target in process.iterfind('targets/targetValue'):
                new_t = rba_xml.TargetSpecies(target.get('species'))
                new_t.value = self._parse_value(
                    target, target.get('species') + '_concentration'
                )
                if target.get('degradation') == '1':
                    new_group.targets.degradation_fluxes.append(new_t)
                elif target.get('dilution_compensation') == '0':
                    new_group.targets.production_fluxes.append(new_t)
                else:
                    new_group.targets.concentrations.append(new_t)
            for target in p.iterfind('targets/targetReaction'):
                new_t = rba_xml.TargetReaction(target.get('reaction'))
                new_t.value = self._parse_value(
                    target, target.get('reaction') + '_flux'
                    )
                new_group.targets.reaction_fluxes.append(new_t)
            if not new_group.is_empty():
                self.model.targets.target_groups.append(new_group)

    def adapt_metabolite_names(self):
        for m in self.model.metabolism.species:
            m.id = metabolite_name(m.id)
        for r in self.model.metabolism.reactions:
            for sr in chain(r.reactants, r.products):
                sr.species = metabolite_name(sr.species)
        for e in self.model.enzymes.enzymes:
            for sr in e.machinery_composition.reactants:
                sr.species = metabolite_name(sr.species)
        for p in self.model.processes.processes:
            mc = p.machinery.machinery_composition
            for sr in chain(mc.reactants, mc.products):
                sr.species = metabolite_name(sr.species)
            for t in chain(p.targets.concentrations,
                           p.targets.production_fluxes,
                           p.targets.degradation_fluxes):
                t.species = metabolite_name(t.species)
        for m in self.model.processes.component_maps:
            for s in chain(m.constant_cost.reactants,
                           m.constant_cost.products):
                s.species = metabolite_name(s.species)
            for c in m.costs:
                for s in chain(c.reactants, c.products):
                    s.species = metabolite_name(s.species)
        for fn in self.model.parameters.functions:
            if fn.variable:
                fn.variable = metabolite_name(fn.variable)

    def add_matp_process():
        result = rba_xml.Process('P_maintenance_atp', 'Maintenance ATP')
        self.model.processes.process.append(result)
        target = rba_xml.TargetReaction('Eatpm')
        target.lower_bound = 'maintenanceATP'
        result.targets.reaction_fluxes.append(target)


class EnzymeConverter(object):
    def __init__(self, xml_root, model, medium):
        self._root = xml_root
        self.model = model
        self._read_enzymes()
        self._read_efficiencies(medium)
        self._convert_metabolite_names()

    def _read_enzymes(self):
        enzymes = self.model.enzymes.enzymes
        for e in self._root.iterfind('listOfEnzymes/enzyme'):
            id_ = e.get('id')
            if self._transport_fns(e):
                forward = self._transporter(id_)
            else:
                forward = self._base(id_)
            new_e = rba_xml.Enzyme(self._name(id_), id_,
                                   forward, self._backward(id_),
                                   rba_xml.is_true(e.get('zero_cost')))
            enzymes.append(new_e)
            n = rba_xml.get_unique_child(e, 'machineryComposition', False)
            if n is not None:
                new_e.machinery_composition \
                    = rba_xml.MachineryComposition.from_xml_node(n)

    def _is_transporter(self, enzyme):
        return enzyme.find('enzymaticActivity/transporterEfficiency/function')

    def _transport(self, reaction_id):
        return reaction_id + '_transport_efficiency'

    def _base(self, reaction_id):
        return reaction_id + '_base_efficiency'

    def _name(self, reaction_id):
        return reaction_id + '_enzyme'

    def _read_efficiencies(medium_id='medium_2'):
        root = self._xml_root('enzymes.xml')
        fn_type = root.find('listOfEfficiencyFunctions/function[@id={}]'
                            .format(medium_id)).get('type')
        for e in root.iterfind('listOfEnzymes/enzyme'):
            id_ = e.get('id')
            base_fn = rba_xml.Function(self._base(id_), fn_type)
            base_fn.parameters = rba_xml.listOfParameters.from_xml_node(
                e.find('enzymaticActivity/enzymeEfficiency[@function={}]'
                       '/listOfParameters'.format(medium_id))
                )
            self.model.parameters.functions.append(base_fn)
            transport = self._transport_fns(e)
            for fn in transport:
                self.model.parameters.functions.append(fn)
            if transport:
                new_agg = rba_xml.Aggregate(self._transporter(id_),
                                            'multiplication')
                for fn in transport:
                    new_agg.function_references.append(
                        rba_xml.FunctionReference(fn.id)
                    )
                self.model.parameters.aggregates.append(new_agg)

    def _transport_fns(self, enzyme):
        result = []
        for i, fn in enumerate(enzyme.iterfind(
                'enzymaticActivity/transporterEfficiency/function')):
            id_ = '{}_transport_factor_{}'.format(enzyme.get('id'), i)
            new_fn = rba_xml.Function(id_, fn.get('type'),
                                      variable=fn.get('variable'))
            new_fn.parameters = rba_xml.listOfParameters.from_xml_node(fn)
            result.append(new_fn)
        return result

    def _convert_metabolite_names(self):
        for e in self.model.enzymes.enzymes:
            for sr in e.machinery_composition.reactants:
                if sr.species == 'm_siroheme':
                    sr.species = 'm_sheme'
                elif sr.species == 'm_biotin':
                    sr.species = 'm_bio'
                elif sr.species == 'm_b6':
                    sr.species = 'm_py5p'


def metabolite_name(old_name):
    """
    Update name to match standards.
    :param old_name: original metabolite name.
    """
    # if species is not a metabolite, ignore it
    if not old_name.startswith('m_'):
        return old_name
    # 1. capitalize the 'm_'
    # 2. add suffix '_c' or change '_xt' to '_e'
    if old_name.endswith('_xt'):
        return old_name[:-3].capitalize() + '_e'
    else:
        return old_name.capitalize() + '_c'


if __name__ == "__main__":
    main()
