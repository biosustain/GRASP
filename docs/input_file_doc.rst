Input excel file specificiation
=============================================

The input file that describes the kinetic model is an `.xlsx` file. 
It has the following sheets:
 
 - `general`: where general information such as model name, parallelization, etc. are defined;
 - `stoic `: where the model's stoichiometric matrix is defined;
 - `mets`: where all the metabolites are specified;
 - `rxns`: where all reactions are specified;
 - `splitRatios`: [TODO Pedro];
 - `poolConst`: [TODO Pedro];
 - `thermo_ineq_constraints`: [TODO Pedro];
 - `thermoRxns`: where the standard Gibbs energies are specified for each reaction in the model (in kJ/mol);
 - `thermoMets`: where the metabolite concentration used to calculate reaction Gibbs energies are specified (in mol/L);
 - `measRates`: where known reactions fluxes are specified (in micromol/gCDW/h);
 - `protData`: where scaled proteomics data is specified;
 - `metsData`: where scaled metabolomics data is specified;
 - `kinetics1`: where enzyme mechanisms and regulation is specified.


Note that the order of the reactions and metabolites in all sheets must be the same.


General
--------------------------------

In the general sheet the second column must be filled by the user with:

 - model name
 - sampling mode: which can be either ORACLE or rejection. ORACLE mode is used to generate models for only one data point, while rejection mode is used to generate models that can explain multiple data points [TODO Pedro plz check];
 - NLP solver: only relevant if the sampling mode is rejection;
 - number of experimental conditions: basically how many data points there are;
 - number of model structures: how many different model structures should be considered, this refers only to regulation though, the reactions and metabolites must be the same in all structures. To specify multiple structures [... TODO Pedro];
 - number of particles: how many models should be generated;
 - parallel mode: whether or not models should be sampled in parallel;
 - number of cores: if parallel mode is ON (set to `1`), how many cores should be used;
 - percentile of alive particles for SMC [TODO Pedro, will this be relevant?];
 - compute robust fluxes: if the fluxes are known only for some reactions in the model, the unknown fluxes can be calculated by setting this option to `1`;
 - compute thermodynamics: if set to 1 the reaction Gibbs free energies will be calculated from the standard Gibbs energies specified in thermoRxns and the metabolite concentrations specified in thermoMets;
 - initial tolerance:  [TODO Pedro];
 - final tolerance:  [TODO Pedro];


stoic
-----------------------------------

In this sheet the reaction stoichiometry is specified, reactions in the rows and metabolites in the columns (it is basically a transposed stoichiometric matrix).


mets
-----------------------------------

This sheet must have 6 columns:

 - ID: the metabolite ids, to be used in the other sheets;
 - metabolite names: the metabolite names;
 - balanced?: whether or not the metabolite is balanced, i.e. it is consumed and produced in the same amount, it's concentration doesn't change;
 - active?: whether or not the metabolite will be modeled;
 - fixed?: whether or not the metabolite should be treated as a constant;
 - measured?: whether or not the metabolite concentration has been measured experimentally;


rxns
-----------------------------------

This sheet must have 5 columns:
 
 - ID: the reaction ids, to be used in the other sheets;
 - reaction names: the reaction names;
 - transportRxn?: whether or not the reaction is a transport reaction;
 - modelled?: whether or not the reaction is part of the model
 - isoenzymes: if there are isoenzymes they must be specified in this column, e.g. if PFK1 and PFK2 are isoenzymes, in the isoenzymes column you should write PFK in PFK1 and PFK2 rows. This is important when the flux through the whole reaction is known but not how much is catalyzed by each isoenzyme individually. By specifying the isoenzymes in this sheet the fraction of flux catalyzed by each isoenzyme individually is sampled randomly for each model.


splitRatios
-----------------------------------

[TODO Pedro]


poolConst
-----------------------------------

[TODO Pedro]


thermo_ineq_constraints
-----------------------------------

[TODO Pedro]



thermoRxns
-----------------------------------

This sheet must have 3 columns:
 
 - rxn: the reaction ids;
 - ∆Gr'_min (kJ/mol): the minimum standard Gibbs energy for each reaction in kJ/mol, tipically calculated as `mean - standard_deviation`
 - ∆Gr'_max (kJ/mol): the maximum standard Gibbs energy for each reaction in kJ/mol, tipically calculated as `mean + standard_deviation`

The standard Gibbs energies can be obtained from `eQuilibrator <http://equilibrator.weizmann.ac.il/>`_.


thermoMets
-----------------------------------

This sheet must have 3 columns:
 
 - mets: the metabolite ids;
 - min (M): the minimum experimental metabolite concentrations in mol/L, typically calculated as `mean - standard_deviation`;
 - max (M): maximum experimental metabolite concentrations in mol/L, typically calculated as `mean + standard_deviation`.

These concentrations are used, together with the standard Gibbs energies in thermoRxns, to calculate each reaction's Gibbs free energy.


measRates
-----------------------------------

This sheet has 3 columns:
 - Fluxes (umol/gdcw/h): 


mets
-----------------------------------
mets
-----------------------------------