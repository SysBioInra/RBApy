
Common
======

Document the code (Stephan)
---------------------------

Write unit/integration tests (Stephan)
--------------------------------------

IDEA: use same composition for every enzyme (Stephan)
-----------------------------------------------------
This would reduce the number of parameters, thus enabling higher quantitative
stability. The idea is that the enzyme "efficiency" would become a compound
parameter reflecting both the weight of the enzyme and its catalytic activity.
Biggest shortcoming would be the necessity for cofactors, but they could be
included in some form. Even if their stoichiometry is messed up it does not
matter that much. Essential is rather can we produce them in the current medium.

The pipeline could propose two branches: one generating the model with full
parameters, one generating the simpler model where the concentrations of all
enzymes are the same.

Automatic rule to remove multiple enzymes catalyzing same reaction? (Stephan)
-----------------------------------------------------------------------------
Duplicating reactions is really annoying and increasing matrix size a lot.
Is there a way to simplify the problem mathematically or can we find a rule
to eliminate some of the enzymes (maybe depending on size and cofactors)?

Data
====

Stoichiometric coefficients in subtilis are different from original model (Stephan)
--------------------------------------
Growth rate loss is 0.66 -> 0.33. We need to renormalize kE by a factor that
reflects new vs old stoichiometry.
kE was originally computed by counting monomers. When available, it was 
multiplied by np, the number of monomers composing enzyme. We should use a
similar procedure.

Transporters vary considerably in subtilis (Stephan)
-------------------------------------------------
examples:
Tnh3b: nh3 -> nh2,hex, 14% variation in Km
Tno2b: . -> hex, 7%
Tqox, Tcta: o2 -> ., 8%
new transporters: Tco2, Ethis, EacpA


Metabolic model of Arabidopsis is incomplete (Stephan)
------------------------------------------------------
It is impossible to synthesize NAD for example.

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
has been implemented (division by an arbitrary number), it might be good to 
go through it again.

Include information retrieved automatically in helper files (Stephan)
--------------------------------------------------------------------
See if and how we should include this information for each helper file.

Readme unclear about how fasta files are working (Caroline)
-----------------------------------------------------------

Protein cofactors like nad/atp/nadp/coA should be removed (Caroline)
--------------------------------------------------------------------
They are not permanent cofactors, they are just used during reaction catalysis.

Remove '_1' for reactions that were not duplicated (Caroline)
-------------------------------------------------------------

Add a fasta file for proteins to complement uniprot (Caroline)
--------------------------------------------------------------

RBA algorithm
=============

Do not allow 0 concentration substrates to enter the cell (Stephan)
-------------------------------------------------------------------
Historical example is 0 concentration trehalose entering periplasm (because 
import fluxes into periplasm do not depend on substrate concentration). 
Trehalose is then transformed to glucose (assumed to be at nonzero 
concentration) and imported into cytoplasm. We should somehow forbid theses
substrates to enter the cell.

Automatic/intelligent preconditionning? (Stephan)
---------------------------
Scaling is a fixed coefficient. If a user uses different units for enzyme
concentrations, it will not work properly. Maybe we should develop a procedure
that automatically scales colums?

Actually, the current scaling seems to be very unstable, it would be urgent
to think about it more thoroughly. Just removing a useless column or sligthly
changing the scaling coefficient can change run times a lot.

Matrix sparsity structure is stable except for mu = 0. We could
also swap rows and columns to improve sparsity as we would only have to
compute row and column order once at the very beginning.

Question: do numerical instabilities come from column scaling, sparsity
structure or some interaction between the two?

For the moment, I put CPLEX's scaling to agressive. This seems sufficient for 
now.

Remove useless columns and rows? (Stephan)
-------------------------------
Removing zero cost enzymes worked pretty well but removing empty reactions
leads to higher solve times. How far should we go. Here is the list of
things that could be removed:
 - useless reactions: reactions that involve only boundary metabolites. It
 seems to make solving worse, so maybe avoid it until we understant why.
 - useless enzyme backward constraints: if an enzyme catalyzes an irreversible
 reaction it does not need to define a backward constraint, it will already be
 included in the lower bound.
 - process machinery columns: only include processes that have a machinery.
However, as long as run times are unstable because of matrix conditionning, it
does not seem to make sense to worry too much about those issues.

Change enzyme-reaction constraints (Stephan)
----------------------------------
See working note for further detail.

Standardize output (Caroline)
----------------------------
As a first attempt, we should at least separate fluxes, enzyme concentrations,
machinery concentrations, etc.
