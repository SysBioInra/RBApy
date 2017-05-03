import os.path
import sys
sys.path.append(os.path.join(sys.path[0],'..'))
import rba
from rba import rba_xml

def flagella_activation():
    process = rba_xml.Process('P_FLAGELLA', 'Flagella activation')
    target = rba_xml.TargetReaction('Th_1')
    target.function_references.append('number_flagella')
    target.function_references.append('flagella_speed')
    target.function_references.append('flagella_h_consumption')
    process.targets.reaction_fluxes.append(target)
    return process

def read_activities(f):
    fns = []
    activity = {}
    for line in f:
        token = line.rstrip('\n').split('\t')
        if token[0] == 'function':
            # format 'function id type'
            fns.append({'id': token[1], 'type': token[2]})
        else:
            # format 'reaction_id fn_id param1_id param1_value ... paramn_value
            reaction = token[0]
            fn = token[1]
            params = iter(token[2:])
            if not(activity.has_key(reaction)): activity[reaction] = {}
            activity[reaction][fn] \
                = {id_: value for id_, value in zip(params,params)}
    return (fns, activity)

def old_name(new_name):
    if new_name == 'R_maintenance_atp':
        return 'Eatpm'
    else:
        return new_name.rsplit('_',1)[0]

if __name__ == "__main__":
    ## generate base XML files with pipeline
    subtilis = rba.PreRba('subtilis/data/params.in')

    ## add flagella constraint
    subtilis.model.processes.processes.append(flagella_activation())

    ## add enzymatic activities
    enzymes = subtilis.model.enzymes
    with open('subtilis/old_data/catalytic_activity.csv', 'r') as f:
        fns, activity = read_activities(f)
    for fn in fns:
        enzymes.efficiency_functions \
               .append(rba_xml.Function(fn['id'], fn['type']))
    for enzyme in enzymes.enzymes:
        reaction = old_name(enzyme.enzymatic_activity.reaction)
        for fn_id in activity[reaction]:
            params = activity[reaction][fn_id]
            new_eff = rba_xml.EnzymeEfficiency(fn_id, params)
            enzyme.enzymatic_activity.enzyme_efficiencies.append(new_eff)

    ## add zero_cost flags
    with open('subtilis/old_data/zero_cost.csv', 'r') as f:
        zero_cost = []
        for line in f:
            zero_cost += line.rstrip('\n').split('\t')
    for enzyme in enzymes.enzymes:
        id_ = old_name(enzyme.id)
        if id_ in zero_cost:
            enzyme.zero_cost = True
            zero_cost.remove(id_)

    ## set medium to original medium
    with open('subtilis/old_data/medium.csv', 'r') as f:
        old_medium = {}
        for line in f:
            met, conc = line.rstrip('\n').split('\t')
            old_medium[met[:-2]] = conc
    subtilis.model.medium = dict.fromkeys(subtilis.model.medium, 0)
    for met, conc in old_medium.iteritems():
        subtilis.model.medium[met] = conc
            
    ## write xml files
    subtilis.model.write_files()
