import sys

from lxml import etree

def dump_all_activities():
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
                f.write('\t'.join([reaction, fn] + params_as_string(efficiency))
                        + '\n')

                
def dump_activities(medium):
    enzymes = etree.ElementTree(file='enzymes.xml')
    with open('catalytic_activity_' + medium + '.csv', 'w') as f:
        xpath = 'listOfEnzymes/enzyme/enzymaticActivity'
        for activity in enzymes.iterfind(xpath):
            reaction = activity.get('reaction')
            for efficiency in activity.iterfind('enzymeEfficiency'):
                if efficiency.get('function') == medium:
                     f.write(reaction + '\t' + params_as_string(efficiency)
                             + '\n')

                     
def params_as_string(efficiency_node):
    params = []
    for p in efficiency_node.iterfind('listOfParameters/parameter'):
        params += [p.get('id'), p.get('value')]
    return '\t'.join(params)
    
    
if __name__ == "__main__":
    if len(sys.argv) == 1:
        dump_all_activities()
    else:
        dump_activities(sys.argv[1])
