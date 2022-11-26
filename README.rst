RBApy
==============================

RBApy is a open-source Python package for the automated generation of bacterial Resource Balance Analysis (RBA) models (https://rba.inrae.fr).
Existing RBA models for *Bacillus subtilis 168* (wild type), *Escherichia coli K-12* (wild type) and CO2-fixing *Escherichia coli K-12* (an engineered strain) can be found here: https://github.com/SysBioInra/Bacterial-RBA-models.

For a complete documentation on RBApy installation and usage, please visit the website:
https://sysbioinra.github.io/RBApy/


Installation
-------------

RBApy requires one of the linear programming solvers `IBM CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_, `GLPK <https://www.gnu.org/software/glpk/>`_, or `Gurobi <https://www.gurobi.com/products/gurobi-optimizer/>`_. Note, while GLPK is capable of executing the example models in the tutorial, in our experience, GLPK is prohibitively slow for real RBA models. IBM and Gurobi both provide free licenses for academic research.

1. Optionally, install the CPLEX linear programming solver.

2. Install this package from PyPI:
    ```
    pip install rbapy
    ```

    Optionally, also install CPLEX by installing RBApy with the ``cplex`` option. Note, this requires a CPLEX license.:

    ```
    pip install rbapy[cplex]
    ```

    Optionally, also install GLPK by installing RBApy with the ``swiglpk`` option. Note, this requires ``libglpk-dev``.:
    
    ```
    pip install rbapy[swiglpk]
    ```

    Optionally, also install Gurobi by installing RBApy with the ``gurobi`` option. Note, this requires a Gurobi license.:

    ```
    pip install rbapy[gurobi]
    ```

More information about how to install RBApy is available at https://sysbioinra.github.io/RBApy/installation.html.

Remark on usage
---------------

When using RBApy in own Python code, it must be imported as: rba

    ```
    import rba
    ```
    
Running
-------

For more detailed instructions on usage, please refer to the RBApy website https://sysbioinra.github.io/RBApy/.

Put the SBML file containing the Genome Size Metabolic Model (GSMM) of your organism of interest in the `input` directory and fill in the
parameter file `input/params.in`. Then open a console at the root
of the repository and run:


```
python generate_rba_model.py input/params.in
```

or, more generally:


```
python generate_rba_model.py path/to/params.in
```

The script will generate several files used as an input for the RBA solver.
By default, they will be written to the `output` directory.

Warning: for the first run, the script will download and parse Uniprot data
as best it can. Unfortunately, numerous values cannot be parsed properly and
are replaced with default values. The script will generate many helper files
to replace these default values with hand-curated values. You should fill in
these helper files and rerun the script to obtain more relevant output
(see instructions below).

If rbapy has been installed properly, there also exists command-line interface called by:


```
generate-rba-model path/to/params.in
```

Once the RBA model has been generated, you can solve it using:


```
python solve_rba_model.py path/to/model
```

where the path points to the directory containing the XML files defining
the RBA model.

If rbapy was installed properly, there also exists command-line interface called by:


```
solve-rba-model path/to/model
```


Running RBApy with Gurobi
^^^^^^^^^^^^^^^^^^^^^^^^^

To use RBApy with Gurobi, either:

* Save your license to your home directory (``~/gurobi.lic``) or to the appropriate location for your OS (e.g., ``/opt/gurobi/gurobi.lic`` for Linux).
* Encode your license variables (e.g., ``WLSACCESSID``, ``WLSSECRET``, ``LICENSEID``) into environment variables with the prefix ``GRB_`` (e.g., ``GRB_WLSACCESSID``, ``GRB_WLSSECRET``, ``GRB_LICENSEID``).


SBML file requirements
----------------------

The SBML file must be a valid SBML file, with gene-reaction associations.
RBApy assumes that the boolean relation is always “or”s of “and”s, e.g. (g1 and g2) or (g3 and g4)
Moreover, the words  “or” and “and” must be written in lowercase letters.
Empty fields in Gene-association will be interpreted as a diffusion reaction.

Authors
-------

Fischer S. , Bulovic A. , Goelzer A. , Bodeit, O., Dinh M.


Citation
---------------

If you use RBApy for scientific publications, please cite Bulović A., Fischer S., Dinh M., Golib F., Liebermeister W., Poirier C., Tournier L., Klipp E., Fromion V., Goelzer A. (2019), "Automated generation of bacterial resource allocation models", Metabolic Engineering 55 (2019) 12–22. `https://doi.org/10.1016/j.ymben.2019.06.001 <https://doi.org/10.1016/j.ymben.2019.06.001>`_.


License
-------

Copyright (c) 2018 INRAE - MaIAGE - France.

RBApy is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RBApy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RBApy.  If not, see <https://www.gnu.org/licenses/>


