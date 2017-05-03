
import os.path
import sys
sys.path.append(os.path.join(sys.path[0],'..'))

import rba
import numpy
from scipy.sparse import coo_matrix

class Reference(object):
    def __init__(self, input_dir):
        ## read data
        self._input_dir = input_dir
        self.A = numpy.loadtxt(self._input('A.txt'), delimiter=',')
        self.Aeq = numpy.loadtxt(self._input('Aeq.txt'), delimiter=',')
        self.b = numpy.loadtxt(self._input('b.txt'), delimiter=',')
        self.beq = numpy.loadtxt(self._input('beq.txt'), delimiter=',')
        self.UB = numpy.loadtxt(self._input('UB.txt'), delimiter=',')
        self.LB = numpy.loadtxt(self._input('LB.txt'), delimiter=',')
        info = self._read_info()
        self.mu = float(info['mu'][0])
        self.metabolites = info['metabolites']
        self.reactions = info['reactions']
        
        ## post processing
        nb_R = len(self.reactions)
        nb_M = len(self.metabolites)
        M_range = {}
        prev = 0
        for name, size in zip(info['metabolite_sets'],
                              info['metabolite_set_size']):
            M_range[name] = prev + numpy.arange(int(size))
            prev += int(size)
        # change signs
        sign_change = 2+numpy.concatenate([M_range[m] for m in ['Xp', 'Xpm', 'Xm', 'Xc']])
        self.Aeq[sign_change] *= -1
        self.beq[sign_change] *= -1
        # move macromolecules to beq
        for i in M_range['Xc']:
            assert(len(numpy.nonzero(self.Aeq[2+i,:])[0])==1)
            reaction_index = numpy.nonzero(self.Aeq[2+i,:])[0][0]
            self.beq -= self.beq[2+i] * self.Aeq[:,reaction_index]
        # remove useless columns and lines
        R_mask = numpy.ones(nb_R, dtype='bool')
        R_mask[-12:] = False
        M_mask = numpy.ones(nb_M, dtype='bool')
        M_mask[M_range['Xc']] = False
        M_mask[self.metabolites.index('hex')] = False
        M_mask[self.metabolites.index('2me4c')] = False
        # rows and colum indices
        enzyme_cols = numpy.arange(nb_R)[R_mask]
        process_cols = nb_R + numpy.arange(2)
        reaction_cols = 2 + nb_R + numpy.arange(nb_R)[R_mask] 
        density_rows = numpy.arange(2)
        forward_rows = 2 + numpy.arange(nb_R)[R_mask]
        backward_rows = 2 + nb_R + numpy.arange(nb_R)[R_mask]
        process_rows = numpy.arange(2)
        metab_rows = 2 + numpy.arange(nb_M)[M_mask]
        reaction_rows = 2 + nb_M + numpy.arange(1)
        # reorder matrices
        rows = numpy.concatenate([forward_rows, backward_rows, density_rows])
        cols = numpy.concatenate([reaction_cols, process_cols, enzyme_cols])
        self.A = self.A[numpy.ix_(rows,cols)]
        self.b = self.b[rows]
        rows = numpy.concatenate([metab_rows, process_rows, reaction_rows])
        self.Aeq = self.Aeq[numpy.ix_(rows, cols)]
        self.beq = self.beq[rows]
        self.UB = self.UB[cols]
        self.LB = self.LB[cols]

    def _read_info(self):
        info = {}
        with open(self._input('info.txt'), 'r') as f:
            for line in f:
                id_, data = line.rstrip('\n').split('=')
                info[id_.strip()] = data.strip().split('\t')
        return info

    def _input(self, file_name):
        return os.path.join(self._input_dir, file_name)

class Candidate(object):
    def __init__(self, input_dir, medium, mu):
        model = rba.RbaModel.from_xml(input_dir)
        solver = model.solve(medium)
        solver.build_matrices(mu)
        self.reactions = solver._model.reactions
        self.metabolites = solver._model.metabolites
        
        # metadata
        nb_R = len(self.reactions)
        nb_P = len(solver._model.processes.ids)
        nb_M = len(self.metabolites)
        R_mask = numpy.ones(nb_R, dtype='bool')
        R_mask[-3:] = False
        P_mask = numpy.ones(nb_P, dtype='bool')
        P_mask[2:] = False
        forward_rows = numpy.arange(nb_R)[R_mask]
        backward_rows = nb_R + numpy.arange(nb_R)[R_mask]
        density_rows = 2*nb_R + numpy.arange(2)
        reaction_cols = numpy.arange(nb_R)[R_mask]
        enzyme_cols = nb_R + numpy.arange(nb_R)[R_mask]
        process_cols = 2*nb_R + numpy.arange(nb_P)[P_mask]
        metab_rows = numpy.arange(nb_M)
        process_rows = nb_M + numpy.arange(nb_P)[P_mask]
        reaction_rows = nb_M + nb_P + numpy.arange(1)
        
        # reorder matrices
        rows = numpy.concatenate([forward_rows, backward_rows, density_rows])
        cols = numpy.concatenate([reaction_cols, process_cols, enzyme_cols])
        self.A = solver.A.toarray()[numpy.ix_(rows,cols)]
        self.b = solver.b[rows]
        self.ineq_row_names = [solver.ineq_row_names[i] for i in rows]
        rows = numpy.concatenate([metab_rows, process_rows, reaction_rows])
        self.Aeq = solver.Aeq.toarray()[numpy.ix_(rows,cols)]
        self.beq = solver.beq[rows]
        self.eq_row_names = [solver.eq_row_names[i] for i in rows]
        self.col_names = [solver.col_names[i] for i in cols]
    
def compare_matrices(new, ref):
    result = []
    diff = coo_matrix(new-ref)
    for i,j,v in zip(diff.row, diff.col, diff.data):
        if ref[i,j] == 0 or abs(v/ref[i,j]) > 1e-3:
            result.append([i, j, v/ref[i,j], new[i,j], ref[i,j]])
    return result

if __name__ == "__main__":
    ref = Reference('ref_data')
    cand = Candidate('subtilis', 'medium_2', ref.mu)
    output = open('subtilis_test.log','w')
    
    ## compare matrices
    ineq_row_names = cand.ineq_row_names
    eq_row_names = cand.eq_row_names
    col_names = cand.col_names
    diff = compare_matrices(cand.A, ref.A)
    for d in diff:
        row, col = d[0], d[1]
        if (ref.UB[row] == 0) and (d[3] == 0):
            continue
        output.write(' '.join([ineq_row_names[row], col_names[col]]
                              + [str(el) for el in d[2:]]) + '\n')
    diff = compare_matrices(numpy.array([cand.b]), numpy.array([ref.b]))
    for d in diff:
        output.write(ineq_row_names[d[1]] + ' '
                     + ' '.join([str(el) for el in d[2:]]) + '\n')
    diff = compare_matrices(cand.Aeq, ref.Aeq)
    for d in diff:
        output.write(' '.join([eq_row_names[d[0]], col_names[d[1]]]
                              + [str(el) for el in d[2:]]) + '\n')
    diff = compare_matrices(numpy.array([cand.beq]), numpy.array([ref.beq]))
    for d in diff:
        output.write(eq_row_names[d[1]] + ' '
                     + ' '.join([str(el) for el in d[2:]]) + '\n')
