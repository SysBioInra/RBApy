from lxml import etree

if __name__ == "__main__":
    enzymes = etree.ElementTree(file='enzymes.xml')
    with open('catalytic_activity.csv', 'w') as f:
        # dump functions
        xpath = 'listOfEfficiencyFunctions/function'
        for fn in enzymes.iterfind(xpath):
            f.write('\t'.join(['function', fn.get('id'), fn.get('type')])+'\n')
        # dump activities
        xpath = 'listOfEnzymes/enzyme/enzymaticActivity'
        for activity in enzymes.iterfind(xpath):
            reaction = activity.get('reaction')
            for efficiency in activity.iterfind('enzymeEfficiency'):
                fn = efficiency.get('function')
                params = {}
                for p in efficiency.iterfind('listOfParameters/parameter'):
                    params[p.get('id')] = p.get('value')
                f.write('\t'.join([reaction, fn] +
                                  [k + '\t' + v for k,v in params.iteritems()])
                        + '\n')
    
