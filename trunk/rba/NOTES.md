
Removing zero_cost enzymes (04/05/17)
=========================

We would like to remove zero_cost enzymes from the current problem. They yield
zero columns that may make the problem longer to solve. Also we wonder how
we can treat enzymes catalyzing reversible reactions as they yield one
unnecessary constraint (backward flux constraint).


Idea 1: Removing them at the end
--------------------------------
We keep indices of spontaneous reactions and zero_cost enzymes. Before every
optimization, we remove the colums corresponding to these enzymes. However,
we also need to remove all the backward and forward flux rows that correspond
to these enzymes. It is not hard but it seems unnecessary.

Idea 2: Adding constraints on the fly for non-zero enzymes
------------------------------------
When we parse enzymes, we can check sort out all enzymes that have a zero
composition (or zero_cost flag set to true). We only build C\_E, W\_E, 
PC\_E and k\_E for the enzymes that are left. We also build an indicator
matrix that relates an enzyme to the reaction it catalyzes (this matrix
could later be reused to tag reactions catalyzed by the same enzyme). In
this first case, we'd still assume that every enzyme constrains both the
forward and backward reaction flux in a symmetrical manner (import term
aside). Schematically, if 1\_R is the indicator matrix, we'd get something like

$$1_R \nu - import_E k_E E \leq 0$$
$$-1_R \nu - k_E E \leq 0$$

Idea 3: Idea 2 + removing useless constraints due to irreversibility
--------------------------------------------------------------------
In this case, we use the same technique as idea 2, except that each enzyme has
to define a forward and a backward constraint (if necessary). There is only one
indicator matrix. The sign of the coefficient in this matrix would reflect
whether we are defining a forward or a backward constraint. 

Advantages are:

 - We would only define one set of constraints instead of two.
 - We can remove all useless backward constraints for enzymes catalyzing
 irreversible reactions.
 
Problems are:
 - Not sure how we define the k\_E part properly, as there is some redundant
 information. We do not want to compute the same efficiency twice. Actually
 we can keep the simple k\_E vector and use indexing to select the same
 coefficient multiple times.
 
Results
-------

I implemented idea 2, leading to a small time improvement (a couple of seconds).
Idea did not work out that well, so I dropped it for the moment. Cplex seems to
handle things decently well. We really need to worry about the stability of the
algorithm.

Removing useless reactions (05/05/17)
==========================

In the E. coli model, there are a lot of reactions that are useless.
I tried removing them: solving time becomes terrible (nearly a 2x increase).
It is not clear to me why the problem would get harder to solve for CPLEX.
I checked that all columns removed where only zeros and the solution found was
the same as the original. Maybe column preprocessing in some form is necessary.
But even if that is true, removing useless reactions mainly improves the
sparsity structure so I'm not sure what puts the solver off.

After further study it appears that removing useless reactions sometimes also
removed useless enzymes. The composition of these enzymes (although useless
in the theoretical problem and solution) seemed to balance the problem from
a numerical point of view... In the end it all breaks down to a problem of
column preprocessing. With a different scaling coefficient (1500 instead
of 1000), everything is fine again. But the system is extremely sensitive to
scaling, which is not good...

