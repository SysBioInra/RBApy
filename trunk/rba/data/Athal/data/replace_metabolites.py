
MISSING = '[MISSING]'

def replace_trnas(lines):
    for l in lines:
        if l[0].startswith('TRNA') and l[2] == MISSING:
            aa = l[0][-3:].capitalize()
            l[2] = 'S_tRNA_40_' + aa + '_41__c'
        elif l[0].endswith('TRNA') and l[2] == MISSING:
            aa = l[0][:3].capitalize()
            l[2] = 'S_L_45__45_tRNA_40_' + aa + '_41__c'

if __name__ == "__main__":
    with open('metabolites.tsv','r') as f:
        lines = [l.rstrip('\n').split('\t') for l in f]
    replace_trnas(lines)
    with open('metabolites.tsv','w') as f:
        f.write('\n'.join(map('\t'.join, lines)))
