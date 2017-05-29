
import copy

import rba

## load data
with open('subtilis/old_data/stoichiometry.csv', 'r') as f:
    data = {}
    for line in f:
        [species, sto] = line.rstrip('\n').split('\t')
        data[species] = float(sto)
model = rba.RbaModel.from_xml('subtilis')

## solve problem
results = []
enzymes = model.enzymes.enzymes
for enzyme in enzymes:
    original_machinery = copy.deepcopy(enzyme.machinery_composition)
    for sr in enzyme.machinery_composition.reactants:
        try:
            sr.stoichiometry = data[sr.species]
        except KeyError:
            print sr.species
    sol = model.solve('medium_2')
    results.append([enzyme.id, sol.mu_opt])
    enzyme.machinery_composition = original_machinery

## print results
results = sorted(results, key=lambda x:-x[1])
with open('subtilis/old_data/stoichiometry_impact.csv', 'w') as f:
    for r in results:
        f.write(r[0] + '\t' + str(r[1]) + '\n')
