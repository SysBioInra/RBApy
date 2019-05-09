

from __future__ import absolute_import, division, print_function

import itertools
import numpy
import scipy.io
import os.path
import collections
import json


class Results(object):

    def __init__(self, model, matrices, solver):
        self._model = model
        self._matrices = matrices
        self._solver = solver
        self.mu_opt = solver.mu_opt
        self.variables = {name: value for name, value in
                          zip(matrices.col_names, solver.X)}
        self.dual_values = {name: value for name, value in
                            zip(matrices.row_names, solver.lambda_)}
        self.reactions = {r.id: reaction_string(r)
                          for r in model.metabolism.reactions}
        self.enzymes = [e.id for e in model.enzymes.enzymes]
        self.processes = [p.id for p in model.processes.processes]
        boundary_metabolites = set(m.id for m in model.metabolism.species
                                   if m.boundary_condition)
        self.boundary_reactions = []
        for r in model.metabolism.reactions:
            if any(m.species in boundary_metabolites
                   for m in itertools.chain(r.reactants, r.products)):
                self.boundary_reactions.append(r.id)

    def reaction_fluxes(self):
        return {i: self.variables[i] for i in self.reactions}

    def enzyme_concentrations(self):
        return {i: self.variables.get(i, 0) for i in self.enzymes}

    def process_machinery_concentrations(self):
        return {i: self.variables[i + '_machinery'] for i in self.processes}

    def sorted_boundary_fluxes(self):
        """
        Return all nonzero import reactions.

        Parameters
        ----------
        solver : RbaSolver
            Solver storing RBA problem solution.

        Returns
        -------
        out : List of (reaction_id, flux) tuples
            List is sorted by reactions with higher flux first.

        """
        return self.sorted_fluxes(self.boundary_reactions)

    def sorted_fluxes(self, ids=None):
        if not ids:
            ids = list(self.reactions.keys())
        result = []
        for reaction in ids:
            flux = self.variables[reaction]
            if flux != 0:
                result.append((self.reactions[reaction], flux))
        result.sort(key=lambda x: abs(x[1]), reverse=True)
        return result

    def saturated_fluxes(self):
        """Return fluxes that are within 90% of upper or lower bound."""
        saturated_lb = {}
        saturated_ub = {}
        for c in self._matrices.reaction_cols:
            name = self._matrices.col_names[c]
            nu = self._solver.X[c]
            ub = self._matrices.UB[c]
            lb = self._matrices.LB[c]
            if lb != 0 and abs((nu-lb)/lb) < 0.1:
                saturated_lb[name] = '{} vs LB {}'.format(nu, lb)
            if ub != 0 and abs((nu-ub)/ub) < 0.1:
                saturated_ub[name] = '{} vs UB {}'.format(nu, ub)
        return saturated_lb, saturated_ub

    def density_status(self, compartment):
        """Return density level for given compartment."""
        row = self._matrices.row_names.index(compartment + '_density')
        max_volume = self._matrices.b[row]
        enzyme_volume = self._matrices.A.toarray()[row, :] * self._solver.X
        occupancy = sum(enzyme_volume)
        enzyme_density = 100 * enzyme_volume / max_volume
        result = [
            (self._matrices.col_names[i], enzyme_volume[i], enzyme_density[i])
            for i in numpy.argsort(enzyme_density)[::-1]
            if enzyme_volume[i] > 0
            ]
        return max_volume, occupancy, result

    def excess_enzymes(self):
        """Find enzymes produced in greater amounts than needed."""
        # find non saturated enzymes
        row_names = self._matrices.row_names
        b_opt = self._matrices.A.dot(self._solver.X)
        result = {}
        for c in self._matrices.enzyme_cols:
            # find capacity rows
            enzyme = self._matrices.col_names[c]
            f = row_names.index(enzyme + '_forward_capacity')
            b = row_names.index(enzyme + '_backward_capacity')
            if (b_opt[f] < -1e-10) and (b_opt[b] < -1e-10):
                result[enzyme] = (b_opt[f], b_opt[b])
        return result

    def write(self, output_dir):
        with open(os.path.join(output_dir, 'reactions.out'), 'w') as f:
            f.write('Reaction\tFlux\n')
            for reaction, flux in self.reaction_fluxes().items():
                f.write('{}\t{}\n'.format(reaction, flux))
        with open(os.path.join(output_dir, 'enzymes.out'), 'w') as f:
            f.write('Enzyme\tConcentration\n')
            for enzyme, conc in self.enzyme_concentrations().items():
                f.write('{}\t{}\n'.format(enzyme, conc))
        with open(os.path.join(output_dir,
                               'process_machineries.out'), 'w') as f:
            f.write('Process\tMachinery Concentration\n')
            for process, conc in \
                    self.process_machinery_concentrations().items():
                f.write('{}\t{}\n'.format(process, conc))


    def write_fluxes(self, output_file, file_type='json', merge_isozyme_reactions=True,
                     only_nonzero=True, remove_prefix=False):
        '''
        Writes the simulation-generated fluxes to a file.

        Parameters
        ----------
        output_file : string
            Path to the output file
        file_type: string
            Either json or csv. In case of csv, values are separated
            by semicolon.
        merge_isozyme_reactions: bool
            If true, fluxes through reactions catalzyed by isoenzymes
            will be added together to a single flux
        only_nonzero: bool
            If true, only non-zero fluxes will be exported
        remove_prefix: bool
            If true, removes the beginning prefix "R_" from the
            reaction's name

        Returns
        -------
        dictionary of flux: values written to the file
        '''
        if file_type not in ('json', 'csv'):
            raise ValueError('file_type can be either json on csv')
        rf = self.reaction_fluxes()
        if remove_prefix:
            rf = {k[2:]: v for k, v in rf.items()}
        if only_nonzero:
            rf = {k: v for k, v in rf.items() if v != 0.}
        if merge_isozyme_reactions:
            rf_merged = collections.defaultdict(float)
            for reac_id, flux_val in rf.items():
                if "duplicate" in reac_id:
                    last_idx = reac_id.index('duplicate') - 1
                    rf_merged[reac_id[:last_idx]] += flux_val
                else:
                    rf_merged[reac_id] += flux_val
            rf = rf_merged
        if file_type == 'json':
            with open(output_file, 'w') as fout:
                fout.write(json.dumps(rf, indent=4))
        elif file_type == 'csv':
            with open(output_file, 'w') as fout:
                fout.write('\n'.join(['{};{}'.format(k, v) for k, v in rf.items()]))
        return rf

    def write_proteins(self, output_file, file_type='csv'):
        '''
        Write simulated protein concentrations to a file.
        output_file: string
            Path to the output file
        file_type: str
            File type, can be csv of json

        Returns
        -------
        Dictionary of protein ID - concentrations written to the file
        '''
        prot_conc = collections.defaultdict(float)
        for enz_id, enz_conc in self.enzyme_concentrations().items():
            if enz_conc == 0:
                continue
            enz = self._model.enzymes.enzymes.get_by_id(enz_id)
            for sp in enz.machinery_composition.reactants:
                prot_conc[sp.species] += sp.stoichiometry * enz_conc
        for mach_id, mach_conc in self.process_machinery_concentrations().items():
            if mach_conc == 0:
                continue
            mach = self._model.processes.processes.get_by_id(mach_id)
            for sp in mach.machinery.machinery_composition.reactants:
                prot_conc[sp.species] += sp.stoichiometry * enz_conc
        if file_type == 'csv':
            with open(output_file, 'w') as fout:
                fout.write('\n'.join(['{}\t{}'.format(k, v) for k, v in prot_conc.items()]))
        elif file_type == 'json':
            with open(output_file, 'w') as fout:
                fout.write(json.dumps(prot_conc, indent=4))
        return prot_conc


    def export_matlab(self, output_dir):
        data = {
            'mu_opt': self._solver.mu_opt, 'X': self._solver.X,
            'A': self._matrices.A, 'b': self._matrices.b,
            'f': self._matrices.f, 'UB': self._matrices.UB,
            'LB': self._matrices.LB, 'ColNames': self._matrices.col_names,
            'Enzymes': self.enzymes, 'Processes': self.processes,
            'RowNames': self._matrices.row_names, 'RowSigns': self._matrices.row_signs,
            'Medium': self._model.medium
        }
        scipy.io.savemat(os.path.join(output_dir, 'rba_matlab.mat'),
                         data, do_compression=True)

    def print_main_transport_reactions(self, number=10):
        print('\nTop {} boundary fluxes:'.format(number))
        for flux in self.sorted_boundary_fluxes()[:number]:
            print('{} {}'.format(*flux))


def reaction_string(reaction):
    reactants = ' + '.join(['{} {}'.format(r.stoichiometry, r.species)
                            for r in reaction.reactants])
    products = ' + '.join(['{} {}'.format(p.stoichiometry, p.species)
                           for p in reaction.products])
    return reaction.id + ': ' + reactants + ' -> ' + products
