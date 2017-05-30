
from lxml import etree
import libsbml
from itertools import chain
import sys

sys.path.append('../..')

import rba
from rba import rba_xml

def metabolite_name(old_name):
    """
    Update name to match standards.
    :param old_name: original metabolite name.
    """
    # if species is not a metabolite, ignore it
    if not old_name.startswith('m_'): return old_name
    # 1. capitalize the 'm_'
    # 2. add suffix '_c' or change '_xt' to '_e'
    if old_name.endswith('_xt'):
        return old_name[:-3].capitalize() + '_e'
    else:
        return old_name.capitalize() + '_c'

def read_metabolism(filename):
    # read sbml
    sbml = libsbml.readSBML(filename)
    metabolism = rba_xml.RbaMetabolism()
    # compartments
    for c in sbml.model.compartments:
        metabolism.compartments.append(rba_xml.Compartment(c.id))
    # metabolites
    for s in sbml.model.species:
        metabolism.species.append(rba_xml.Species(s.id, s.boundary_condition))
    # reactions
    for r in sbml.model.reactions:
        new_reaction = rba_xml.Reaction(r.id, r.reversible)
        metabolism.reactions.append(new_reaction)
        for sr in r.reactants:
            new_reaction.reactants.append \
                (rba_xml.SpeciesReference(sr.species, sr.stoichiometry))
        for sr in r.products:
            new_reaction.products.append \
                (rba_xml.SpeciesReference(sr.species, sr.stoichiometry))
    return metabolism

def read_macromolecules(filename, tag):
    root = etree.ElementTree(file=filename).getroot()
    result = rba_xml.RbaMacromolecules()
    # read components
    comp = rba_xml.get_unique_child(root, 'listOfComponents')
    result.components = rba_xml.ListOfComponents.from_xml_node(comp)
    # read macromolecules
    molecules = root.find('listOfSpecies').findall(tag)
    for m in molecules:
        result.macromolecules.append(rba_xml.Macromolecule.from_xml_node(m))
    return result

def read_enzymes(filename):
    root = etree.ElementTree(file=filename).getroot()
    result = rba_xml.RbaEnzymes()
    # read list of efficiency functions
    eff_fns = rba_xml.get_unique_child(root, 'listOfEfficiencyFunctions')
    result.efficiency_functions \
        = rba_xml.ListOfEfficiencyFunctions.from_xml_node(eff_fns)
    # read enzymes
    enzymes = root.findall('listOfEnzymes/enzyme')
    for e in enzymes:
        new_e = rba_xml.Enzyme(e.get('id'), rba_xml.is_true(e.get('zero_cost')))
        result.enzymes.append(new_e)
        # machinery
        n = rba_xml.get_unique_child(e, 'machineryComposition', False)
        if n is not None:
            new_e.machinery_composition = rba_xml.MachineryComposition.from_xml_node(n)
        # enzyme efficiencies
        effs = e.findall('enzymaticActivity/enzymeEfficiency')
        for eff in effs:
            new_e.enzymatic_activity.enzyme_efficiencies.append \
                (rba_xml.EnzymeEfficiency.from_xml_node(eff))
        # transport efficiency
        n = e.find('enzymaticActivity/transporterEfficiency')
        if n is not None:
            new_e.enzymatic_activity.transporter_efficiency \
                = rba_xml.TransporterEfficiency.from_xml_node(n)
    return result

def read_processes(filename):
    root = etree.ElementTree(file=filename).getroot()
    result = rba_xml.RbaProcesses()
    
    ## processes
    processes = root.find('listOfProcesses')
    for p in processes:
        new_p = rba_xml.Process(p.get('id'), p.get('name'))
        result.processes.append(new_p)
        # machinery
        n = p.find('capacityConstraint')
        if n is not None: new_p.machinery = rba_xml.Machinery.from_xml_node(n)
        # operations
        productions = p.findall('operatingCosts/production')
        for prod in productions:
            new_p.operations.productions.append \
                (rba_xml.Operation.from_xml_node(prod))
        degradations = p.findall('operatingCosts/degradation')
        for deg in degradations:
            new_p.operations.degradations.append \
                (rba_xml.Operation.from_xml_node(deg))
        # targets
        targets = p.findall('targets/targetValue')
        for t in targets:
            new_t = rba_xml.TargetSpecies(t.get('species'))
            # read value and function references
            if t.get('value'): new_t.value = float(t.get('value'))
            fn_refs = t.findall('functionReference')
            for fn_ref in fn_refs:
                new_t.function_references.append(fn_ref.get('function'))
            # detect whether target is concentration or absolute flux.
            if t.get('degradation') == '1':
                new_p.targets.degradation_fluxes.append(new_t)
            elif t.get('dilution_compensation') == '0':
                new_p.targets.production_fluxes.append(new_t)
            else:
                new_p.targets.concentrations.append(new_t)
        targets = p.findall('targets/targetReaction')
        for t in targets:
            new_t = rba_xml.TargetReaction(t.get('reaction'))
            # read value and function references
            if t.get('value'): new_t.value = float(t.get('value'))
            fn_refs = t.findall('functionReference')
            for fn_ref in fn_refs:
                new_t.function_references.append(fn_ref.get('function'))
            new_p.targets.reaction_fluxes.append(new_t)

    ## component maps
    maps = root.findall('listOfComponentMaps/componentMap')
    for m in maps:
        new_m = rba_xml.ComponentMap(m.get('id'))
        result.component_maps.append(new_m)
        n = m.find('constantCost')
        if n is not None:
            new_m.constant_cost = rba_xml.ConstantCost.from_xml_node(n)
        costs = m.findall('cost')
        for c in costs:
            new_m.costs.append(rba_xml.Cost.from_xml_node(c))
    return result

if __name__ == "__main__":
    output = '.'
    input_ = 'old_data/'
    model = rba.RbaModel()
    
    ## convert old files
    # metabolism
    model.metabolism = read_metabolism(input_ + 'metabolism.xml')
    # parameters
    model.parameters = rba_xml.RbaParameters.from_file(input_+'parameters.xml')
    # macromolecules
    model.proteins = read_macromolecules(input_ + 'proteins.xml', 'protein')
    model.rnas = read_macromolecules(input_ + 'rnas.xml', 'rna')
    model.dna = read_macromolecules(input_ + 'dna.xml', 'dna')
    # enzymes
    model.enzymes = read_enzymes(input_ + 'enzymes.xml')
    # change metabolite names where necessary
    for e in model.enzymes.enzymes:
        for sr in e.machinery_composition.reactants:
            if sr.species == 'm_siroheme': sr.species = 'm_sheme'
            elif sr.species == 'm_biotin': sr.species = 'm_bio'
            elif sr.species == 'm_b6': sr.species = 'm_py5p'
    # processes
    model.processes = read_processes(input_ + 'processes.xml')
    # add maintenance atp process
    p_atpm = rba_xml.Process('P_maintenance_atp', 'Maintenance ATP')
    target = rba_xml.TargetReaction('Eatpm')
    target.function_references.append('maintenanceATP')
    p_atpm.targets.lower_bounds.append(target)
    model.processes.processes.append(p_atpm)
    # medium
    model.medium = model.read_medium('medium.tsv')

    ## adapt metabolite names
    # metabolism
    for m in model.metabolism.species:
        m.id = metabolite_name(m.id)
    for r in model.metabolism.reactions:
        for sr in chain(r.reactants, r.products):
            sr.species = metabolite_name(sr.species)
    # enzymes
    for e in model.enzymes.enzymes:
        for sr in e.machinery_composition.reactants:
            sr.species = metabolite_name(sr.species)
        for fn in e.enzymatic_activity.transporter_efficiency:
            fn.variable = metabolite_name(fn.variable)
    # processes
    for p in model.processes.processes:
        mc = p.machinery.machinery_composition
        for sr in chain(mc.reactants, mc.products):
            sr.species = metabolite_name(sr.species)
        for t in chain(p.targets.concentrations, p.targets.production_fluxes,
                       p.targets.degradation_fluxes):
            t.species = metabolite_name(t.species)
    for m in model.processes.component_maps:
        for s in chain(m.constant_cost.reactants, m.constant_cost.products):
            s.species = metabolite_name(s.species)
        for c in m.costs:
            for s in chain(c.reactants, c.products):
                s.species = metabolite_name(s.species)
            
    ## write model to file
    model.write_files(output)
