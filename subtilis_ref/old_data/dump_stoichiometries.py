from lxml import etree

if __name__ == "__main__":
    enzymes = etree.ElementTree(file='enzymes.xml')
    xpath = 'listOfEnzymes/enzyme/machineryComposition/listOfReactants/speciesReference'
    data = {}
    inconsistent = set()
    for protein in enzymes.iterfind(xpath):
        species = protein.get('species')
        if species.startswith('m_'): continue
        sto = protein.get('stoichiometry')
        try:
            if data[species] != sto:
                inconsistent.add(species)
        except KeyError:
            data[species] = sto
    print('Inconsistent data for species ' + ', '.join(inconsistent))
    with open('stoichiometry.csv', 'w') as f:
        f.write('\n'.join(['\t'.join([spe, sto]) for spe, sto in data.iteritems()]))
    
