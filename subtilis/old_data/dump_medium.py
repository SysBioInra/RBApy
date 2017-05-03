import libsbml

if __name__ == "__main__":
    # read sbml
    sbml = libsbml.readSBML('sbml.xml')
    medium = {}
    for m in sbml.model.species:
        if m.compartment == 'e' and m.initial_concentration > 0:
            medium[m.id] = m.initial_concentration
    with open('medium.csv', 'w') as f:
        for id_, value in medium.iteritems():
            f.write(id_ + '\t' + str(value) + '\n')
    
