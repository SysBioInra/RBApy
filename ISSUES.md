
Common
======

Document the code (Stephan)
---------------------------

Write unit/integration tests (Stephan)
--------------------------------------

Data
====

Stoichiometric coefficients in subtilis are different from original model (Stephan)
--------------------------------------
Growth rate loss is 0.66 -> 0.33

PreRBA
======

Accept multiple mRNAs to compute average mRNA (Stephan)
-------------------------------------------------------

Accept fasta file with DNA composition (Stephan)
------------------------------------------------

Create a clear design (Stephan)
-------------------------------

Mostly it would be nice to separate

 - parsing, where parts of XML are created on the fly while the rest is stored
 as temporary data used later.
 - post processing, where the temporary data is used to populate the XML with
 enzymes, processes and so on.
 
If we manage to make these structures clear, we can create an API where the
user will be able to write scripts that can be reused for several organisms.
For example, the temporary data contains the metabolite table, that links all
key metabolites to a constant name. A script using the metabolite table would
be portable, whereas a script using metabolite names for a specific organism
(reusing SBML names) could not necessarily be reused as such for a different
organism.

All average enzymes are the same (Stephan)
------------------------------------------

Compared to Anne's original version, all average proteins have the same
composition. Anne had three different proteins:
 
 - average enzymatic protein in the cytosol.
 - average enzymatic protein in the membrane.
 - average non-enzymatic protein.
 
Is the new version actually fine?

Separate enzyme composition and catalytic functions? (Stephan)
--------------------------------------------------------------
Do we need two different files for clarity ?

Add SBML reference + date of extraction somewhere in the model (Ana)
--------------------------------------------------------------------

Target concentrations for tRNAs are wrong (Stephan)
---------------------------------------------------
Targets were originally expressed as a nucleotide concentration, not a tRNA
concentration. We need to divide them by the length of tRNAs. A quick fix
has been implemented, it might be good to go through it again.

RBA algorithm
=============

Bug: link enzymes with reactions properly (Stephan)
---------------------------------------------------

In all test cases, reactions and enzymes are in the exact same order. Enzymes
should be reordered according to reactions they catalyze.

Do not allow 0 concentration substrates to enter the cell (Stephan)
-------------------------------------------------------------------

Historical example is 0 concentration trehalose entering periplasm (because 
import fluxes into periplasm do not depend on substrate concentration). 
Trehalose is then transformed to glucose (assumed to be at nonzero 
concentration) and imported into cytoplasm. We should somehow forbid theses
substrates to enter the cell.

Apply scaling on subblocks? (Stephan)
-------------------------------------

Instead of applying column scaling before solving, we could apply it on 
subblocks so we do it only once.

Automatic scaling? (Stephan)
---------------------------

Scaling is a fixed coefficient. If a user uses different units for enzyme
concentrations, it will not work properly. Maybe we should develop a procedure
that automatically scales colums?

Remove empty rows and columns (Stephan)
---------------------------------------

The algorithm should detect and remove empty rows and columns. Not sure 
how much time we can win with this as CPLEX probably already does that.
For enzymes, this means that we will lose the nice diagonal structure
of the capacity subblock. I am not sure whether it is easier:

 - to build the whole diagonal subblock, then remove columns.
 - populate a matrix coefficient by coefficient.
 

Idea: fusing identical enzymes (Stephan)
----------------------------------

The idea is to go back to $\nu / k_E \leq E$ instead of $\nu \leq k_EE$. 
Thanks to
this change, we can put all reactions related to a single enzyme (*e.g.*
average enzyme, nonspecific transporter) on the same matrix row. 
The row would then read

$$\sum \nu_i / k_{E_i} \leq \sum E_i = E$$

In my opinion, CPLEX should be able to renormalize the line properly. If not,
multiplying all rows by an average k~E~ should do the trick.

