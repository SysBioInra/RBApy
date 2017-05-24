
import xml.etree.cElementTree as etree
import libsbml
import itertools

def metabolite_name(old_name):
    """
    Update name to match standards.
    :param old_name: original metabolite name.
    """
    # 1. capitalize the 'm_'
    # 2. add suffix '_c' or change '_xt' to '_e'
    if old_name.endswith('_xt'):
        return old_name[:-3].capitalize() + '_e'
    else:
        return old_name.capitalize() + '_c'

def enzyme_note(enzyme_list):
    result = '<body xmlns="http://www.w3.org/1999/xhtml"><p>GENE_ASSOCIATION: ('
    result += ' and '.join(enzyme_list)
    return result + ')</p></body>\n'

if __name__ == "__main__":
    # read sbml
    sbml = libsbml.readSBML('metabolism.xml')
    # update metabolite names
    for s in sbml.model.species:
        s.setId(metabolite_name(s.id))
    for r in sbml.model.reactions:
        for sr in itertools.chain(r.reactants, r.products):
            sr.setSpecies(metabolite_name(sr.species))
    # add notes with enzyme compositions
    enzymes = etree.parse('enzymes.xml')
    enzyme_list = {}
    for e in enzymes.iterfind('listOfEnzymes/enzyme'):
        id_ = e.attrib['id']
        enzyme_list[id_] = [sr.attrib['species'] for sr in e.find('machineryComposition/listOfReactants') if not sr.attrib['species'].startswith('m_')]
    for r in sbml.model.reactions:
        r.setNotes(enzyme_note(enzyme_list[r.id]))
    # write sbml
    libsbml.writeSBML(sbml, 'sbml.xml')
