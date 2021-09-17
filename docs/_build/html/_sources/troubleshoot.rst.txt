Troubleshooting
=======================

There are multiple errors that you can get when running GRASP. |br|
Here we go through some of the most common (or tricky) errors, what they mean and what you can do about them.

.. contents:: :local:
    :depth: 3


'The jacobian matrix is not square. Hint: check if the correct metabolites are set as constants.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error can occur when running the function ``checkStability``, which calculates the eigenvalues of the system's Jacobian and checks if the real part of the eigenvalues is below a certain threshold.

The error happens when the jacobian matrix is not square. This is generally because the user has defined the wrong metabolites as balanced in the ``mets`` sheet of the excel input file.


'Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error occurs when using TMFA to find the bounds for Gibbs energies, reaction fluxes, and metabolite concentrations.

It often happens because the lower or upper bound value in the ``thermoMets`` sheet of the excel input file has been set to 0. If a metabolite has concentration of 0 mol/L, a very small value must be used instead, otherwise we will have numerical issues.


'The initial point obtained from TMFA is not thermodynamically feasible.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use TMFA to find ranges for fluxes, Gibbs energies, and metabolite concentrations compatible with each other. We do so by maximizing and minimizing each variable (reaction flux, Gibbs energy, and metabolite concentration) subject to the thermodynamic constraints set by TMFA.

At the end we take the solution to the minimization/maximization of the last variable and take the mean of the miniminum and maximum solution values for each variable as an initial point for the subsequent sampling of Gibbs energies and reactions fluxes that are thermodynamically feasible.

It shouldn't happen, but it is possible that this initial point is not thermodynamically feasible. A point that is thermodynamically feasible satisfies the TMFA constraints and the sum of Gibbs energies for all reactions involved in a closed loop is zero.


'The flux directions are not consistent. Please make sure that both the lower and upper bound of the flux ranges (fluxMean - 2*fluxStd and fluxMean + 2*fluxStd, respectively) are either positive or negative.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error can happen in two situations:

 - when the input excel file is loaded and the flux lower bounds don't have the same sign, where the bounds are calculated as:

   - lower bound: fluxMean - 2*fluxStd
   - upper bound: fluxMean + 2*fluxStd
 - if the flux lower and upper bounds calculated in ``computeGibbsFreeEnergyRanges`` by using TMFA don't have the same sign, i.e. are not both positive or both negative.

In practice it should only happen in the first situtation.

An easy fix is to make sure that the mean and standard deviation given in ``measRates`` in the excel input file lead to lower and upper bounds that are either both positive or both negative.


'The TMFA problem is infeasible. Verify that the standard Gibbs free energy and metabolite concentration values are valid/correct.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you get this error it means that the given ranges for standard Gibbs energies (in ``thermoRxns``), the ranges for reactions fluxes (given in ``measRates``), and the ranges for metabolite concentrations (given in ``thermoMets``) are not compatible with each other.

Along with this error message you also get a list of the metabolite concentrations as well as the reactions whose standard Gibbs energies are causing the error.

You only need to fix either the metabolite concentrations or the standard Gibbs energies, not both.

By using the jupyter notebook, ``visualize_initial_dgs_and_fluxes`` in ``visualization`` you can get an idea of how off the Gibbs energies and reactions fluxes are.


'The TMFA is infeasible. The MILP model has been written in the examples folder so that you can inspect it.'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes the TMFA problem is infeasible and our diagnostics function cannot find which reaction or metabolite is the problem, in that case, we use ``gurobi_iis`` to find which constraint or variable bounds are causing the issue and add that to the error message. 

We also write the MILP problem into the folder where you're running your script from (maybe the ``examples`` folder) so that you can inspect it and find out what the issue is with the given constraint or variable bounds.

.. |br| raw:: html

      <br>