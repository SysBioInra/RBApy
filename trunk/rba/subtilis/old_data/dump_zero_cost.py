from lxml import etree

if __name__ == "__main__":
    enzymes = etree.ElementTree(file='enzymes.xml')
    xpath = 'listOfEnzymes/enzyme'
    zero_cost = []
    for enzyme in enzymes.iterfind(xpath):
        if enzyme.get('zero_cost') == '1':
            zero_cost.append(enzyme.get('id'))
    with open('zero_cost.csv', 'w') as f:
        f.write('\t'.join(zero_cost))
    
