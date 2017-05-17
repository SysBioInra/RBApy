preRBA: Preprocessing Annotation for Resource Balance Analysis
==============================================================

A workflow that obtains various annotation info
from the `UniProt database <https://www.uniprot.org>`_ for any organism.
This module also exports various files that are used for  
initialization of a RBA problem as proposed by
`MaiAGE <http://maiage.jouy.inra.fr>`_.

Getting Started
---------------

This workflow was designed with Python 2.7. Make sure you have
 - the full preRBA repository.
 - a SBML file for your favorite organism (see requirements below).

Running
-------

Put the SBML file in the `input` directory and fill in the
parameter file `input/params.in`. Then open a console at the root
of the repository and run:

```
python preRBA.py input/params.in
```

or, more generally:

```
python preRBA.py path/to/params.in
```
The script will generate several files used as an input for the RBA solver.
By default, they will be written to the `output` directory.

Warning: for the first run, the script will download and parse Uniprot data
as best it can. Unfortunately, numerous values cannot be parsed properly and
are replaced with default values. The script will generate many helper files
to replace these default values with hand-curated values. You should fill in
these helper files and rerun the script to obtain more relevant output
(see instructions below).

Helper Files
------------

Several helper files are generated by the program. They will all be located
in the input directory *after* the script has been run for the first time.
 - `uniprot.csv`: protein data downloaded from uniprot. Users can modify
   this file (but it is usually not necessary) or replace it with a uniprot
   file of their own (at their own risk if important columns are missing!).
 - `unknown_proteins.tsv`: proteins used to define enzymes in the SBML file
   that could not be retrieved in uniprot. User should replace values tagged
   as `[MISSING]` with a uniprot identifier, an average protein (ids of the
   type `average_protein_xxx` as defined in `proteins.xml`) or leave empty
   for spontaneous reactions.
 - `cofactors.tsv`: for every protein, lists name, chebi identifier and
   stoichiometry of cofactors. Based on the last column (containing uniprot
   annotation), user should check inferred values and fill in fields tagged
   as `[MISSING]`. Important: if annotation says that there is *no* cofactor,
   *do not* remove line from csv file. Ignore name and chebi fields and set
   stoichiometry to 0.
 - `subunits.tsv`: this file is used to retrieve the stoichiometry of a
   protein within its enzymatic complex. For every ambiguous entry, gene
   name and protein name are provided for proper identification. User
   should read annotation field and provide `[MISSING]` stoichiometry.
 - `location_mapping.tsv`: this file is used to link uniprot location ids
   with user-defined compartment ids. For every
   uniprot location, user may fill in a name compartment. These may or
   may not be the compartments defined in the SBML. The idea behind this
   file is that it can be used to fuse compartments or ignore them.
   If the field is left empty (or `[MISSING]`), a compartment named after
   uniprot location is created.
 - `locations.tsv`: this file lists proteins for which uniprot has no
   location data. Protein and gene names are given to help fill in missing
   locations. Note that locations *must* correspond to names defined in
   `location_mapping.tsv`. If fields are left empty or `[MISSING]`,
   proteins are assumed to be located in `Cytoplasm` (or user-defined
   couterpart as defined in `location_mapping.tsv`.

In order for the RBA results to be relevant, please fill in as much
information as possible. Every time you provide new information, please
rerun the pipeline to regenerate RBA input properly.

SBML file requirements
----------------------

TODO

Repository structure and files:
-------------------------------

pipeline
  |
  +------ input (default input directory)
  +------ input/params.in (sample parameter file)
  |
  +------ output (default output directory)
  |
  +------ src (source files)
  |
  |
  + preRBA.py (main script)

Authors
-------

License
-------

Acknowledgments
---------------