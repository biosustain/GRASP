How to build a model
=============================================

To build a model you need to fill in an excel file following a specific format. This is the input for GRASP, which will take it and generate a model ensemble.

The input file contains the following information:

 - reaction stoichiometry;
 - reaction fluxes (in mmol/L/h);
 - Gibbs standard energies (in kJ/mol);
 - metabolite concentrations (in mol/L), used to calculate the reactions free Gibbs energies;
 - enzyme mechanisms;
 - enzyme regulation.

To define an enzyme mechanism, a pattern file needs to be specified, which defines the metabolites binding/release order. Instructions on how to specify such file can be found below in :ref:`pattern_file_specification`.

Setting up the excel input file is a fairly tedious and error prone process, thus we advise you to use the python package  `set_up_grasp_models <https://github.com/biosustain/set_up_grasp_models>`_ as a starting point.
This package can also be used to make sure the excel input file is valid.

Below is a description of every sheet and every column that must be filled in the excel input file.
Example input files can also be found in the folder ``io/input``.

GRASP can be used together with a rejection sampler (following the ABC approach). In that case the excel input file will include more data.



Input excel file description
--------------------------------------------

The input file that describes the kinetic model is an ``.xlsx`` file.

It has the following sheets:

 - general_: where general information such as model name, number of models in the ensemble, parallelization, etc. are defined;
 - stoic_: where the model's (transposed) stoichiometric matrix is defined;
 - mets_: where all the metabolites are specified;
 - rxns_: where all reactions are specified;
 - poolConst_: where conservation relations can be defined for groups of pool metabolites (e.g., NAD/NADH) for each experimental condition;
 - thermo_ineq_constraints_: where ratio inequality constraints on metabolite concentrations can be defined;
 - thermoRxns_: where the standard Gibbs energies are specified for each reaction in the model (in kJ/mol);
 - thermoMets_: where the metabolite concentrations used to calculate reaction Gibbs energies are specified (in mol/L);
 - measRates_: where known reactions fluxes are specified (in mmol/L/h);
 - protData_: where scaled proteomics data is specified - this is only relevant if you have more than one data point, i.e. if you are using the ABC approach;
 - metsData_: where scaled metabolomics data is specified - this is only relevant if you have more than one data point, i.e. if you are using the ABC approach;
 - kinetics1_: where enzyme mechanisms and regulation is specified.


Note that the order of the reactions and metabolites in all sheets must be the same.


general
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the general sheet the second column must be filled by the user with:

 - ``model name``: the name of the model;
 - ``sampling mode``: which can be either ``GRASP`` or ``rejection``. GRASP mode is used to generate models for only one data point, while rejection mode is used to generate models that can explain multiple data points;
 - ``NLP solver``: only relevant if the sampling mode is rejection. It can be either ``fmincon`` or ``nlopt``, the latter is more efficient but you need to install it from `here <https://nlopt.readthedocs.io/en/latest/>`_;
 - ``LP solver``: the solver for linear programs, it can be either ``linprog``, the LP solver provided by Matlab, or ``gurobi`` if you have installed it;
 - ``prior distribution for fluxes``: it can be either ``normal`` or ``uniform`` and it is the distribution used when sampling flux values with the Hit-and-Run algorithm;
 - ``prior distribution for thermodynamic quantities``: it can be either ``normal`` or ``uniform`` and it is the distribution used when sampling standard Gibbs free energy and metabolite concentration values with the Hit-and-Run algorithm;
 - ``number of experimental conditions``: how many data points there are besides the reference state. If you are only using GRASP this should be set to ``0``, because there is only the reference state;
 - ``number of model structures``: how many different model structures should be considered, this refers only to regulation though, the reactions and metabolites must be the same in all structures;
 - ``number of models``: how many models should be generated;
 - ``parallel mode``: whether (``1``) or not (``0``) models should be sampled in parallel;
 - ``number of cores``: if parallel mode is ON (set to ``1``), how many cores should be used;
 - ``compute robust fluxes``: if the fluxes are known only for some reactions in the model, the unknown fluxes can be calculated by setting this option to ``1``. For more info on how fluxes are calculated see :ref:`compute_robust_fluxes`;
 - ``rejection sampler tolerance``: this is the maximum deviation of the sampled model to the experimental data allowed in the rejection sampler. Only relevant when using the ABC approach for multiple data points.


stoic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this sheet the reaction stoichiometry is specified, reactions in the rows and metabolites in the columns (it is basically a transposed stoichiometric matrix).


mets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet must have 3 columns:

 - ``metabolite ID``: the metabolite ids, to be used in the other sheets;
 - ``metabolite name``: the metabolite names;
 - ``balanced?``: whether or not the metabolite is balanced, i.e.  whether or not it is consumed and produced in the same amount. If a metabolite is defined as not balanced, then its concentration will be considered to be constant;


rxns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet must have 4 columns:
 
 - ``reaction ID``: the reaction ids, to be used in the other sheets;
 - ``reaction name``: the reaction names;
 - ``transport reaction?``: whether or not the reaction is a transport reaction. This is used only to exclude transport reactions from the TMFA problem;
 - ``isoenzymes``: if there are isoenzymes they must be specified in this column, e.g. if PFK1 and PFK2 are isoenzymes, in the isoenzymes column you should write PFK in PFK1 and PFK2 rows. This is important when the flux through the whole reaction is known but not how much is catalyzed by each isoenzyme individually. By specifying the isoenzymes in this sheet the fraction of flux catalyzed by each isoenzyme individually is sampled randomly for each model.


poolConst
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet is optional and it is only required if pool concentration data is to be used (rejection mode). If the latter is the case, the sheet must have the first column indicating the names of all the active species (metabolites and enzymes, in this order) in the model including an additional row (the ID of the last row is not important, but should be filled). After the first column is defined, there should be one column per experimental condition (excluding the reference condition). The headers of these columns are irrelevant, but should be filled. In each of these columns, there should be a ‘1’ in the position of each conserved species and a ‘0’ otherwise. Finally, the last value in these columns should be filled with the experimental relative pool concentration value for each condition. Currently, GRASP supports only the addition of one of these constraints per experimental condition.

See tutorial_02_adenine_cycle under `io/input` for an example.


thermo_ineq_constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet is optional and it is only required if concentration ratio inequalities are to be used (rejection mode). In log space, these constraints translate into linear relationships that can be included in GRASP. If this information is available, the sheet must have the first column indicating the names of all the considered metabolites in the model including an additional row (the ID of the last row is not important, but should be filled). After the first column is defined, there should be one column per ratio condition. The headers of these columns are irrelevant, but should be filled. In each of these columns, there should be a ‘1’ (metabolite in the numerator) or ‘-1’ (metabolite in the denominator) in the position of the metabolites whose ratios are known. The rest of entries should be filled with ‘0’. Finally, the last value in these columns should be filled with the right-hand side value of the inequality with ‘<’. GRASP treats these constraints globally, which means that they will included in the fitting of all experimental conditions.


thermoRxns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet must have 3 columns:
 
 - ``reaction ID``: the reaction ids;
 - ``∆Gr'_min (kJ/mol)``: the minimum standard Gibbs energy for each reaction in kJ/mol, tipically calculated as ``mean - 2*standard_deviation``
 - ``∆Gr'_max (kJ/mol)``: the maximum standard Gibbs energy for each reaction in kJ/mol, tipically calculated as ``mean + 2*standard_deviation``

The standard Gibbs energies can be obtained from `eQuilibrator <http://equilibrator.weizmann.ac.il/>`_.


thermoMets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet must have 3 columns:
 
 - ``metabolite ID``: the metabolite ids;
 - ``min (M)``: the minimum experimental metabolite concentrations in mol/L, typically calculated as ``mean - 2*standard_deviation``;
 - ``max (M)``: the maximum experimental metabolite concentrations in mol/L, typically calculated as ``mean + 2*standard_deviation``.

These concentrations are used, together with the standard Gibbs energies in ``thermoRxns``, to calculate each reaction's Gibbs free energy.


measRates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet has 3 columns:

 - ``reaction ID``: the reaction ids;
 - ``vref_mean (mmol/L/h)``: the average flux for the reactions whose flux is known. It should be specified in mmol/L/h;
 - ``vref_std (mmol/L/h)``: the standard deviation of the measured flux. It should be specified in mmol/L/h.

Note that here you should only specify fluxes whose values you know and are non-zero. Reactions with zero flux should not be included in the model, as these cannot be parameterized.

Furthermore, the standard deviation should never be zero.


protData
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet is only relevant if there are experimental conditions and the sampling mode is ``rejection``. It has 4 columns:

 - ``reaction / enzyme ID``: the reaction ID
 - ``lower_bound``: the scaled lower bound for the enzyme concentration, typically ``(mean - 2*std) / mean``;
 - ``mean``: the scaled mean value for the enzyme concentration, ``mean/mean``;
 - ``upper_bound``: the scaled upper bound for the enzyme concentration, typically ``(mean + 2*std) / mean``.


metsData
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet is only relevant if there are experimental conditions and the sampling mode is ``rejection``. It sheet has 4 columns:

 - ``metabolite ID``: the metabolite ID
 - ``lower_bound``: the scaled lower bound for the metabolite concentration, typically ``(mean - 2*std) / mean``;
 - ``mean``: the scaled mean value for the metabolite concentration, ``mean/mean``;
 - ``upper_bound``: the scaled upper bound for the metabolite concentration, typically ``(mean - 2*std) / mean``.


kinetics1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sheet has at least 11 columns (more can be added to add notes regarding references and etc.):

 - ``reaction ID``: the reaction IDs
 - ``kinetic mechanism``: the kinetic/enzyme mechanism for the reaction, e.g. ordered Bi Bi. This should be the name of the pattern file in the ``patterns`` folder where the mechanism is specified.
 - ``substrate order``: the binding order for the substrates. Substrates must be separated by a single space.
 - ``product order``: the release order for the products. Products must be separated by a single space.
 - ``promiscuous``: if the enzyme that catalyzes a given reaction is promiscuous this column should include the IDs of all the reactions catalyzed by that enzyme. Reaction IDs must be separated by a single space.
 - ``inhibitors``: metabolite IDs for inhibitors of the given reaction. The binding/release of these inhibitors should be included in the pattern file that describes the enzyme mechanism.
 - ``activators``: metabolite IDs for activators of the given reaction. The binding/release of these activators should be included in the pattern file that describes the enzyme mechanism.
 - ``negative effectors``: metabolite IDs for negative allosteric effectors. These do not need to be included in the reaction's pattern file.
 - ``positive effectors``: metabolite IDs for positive allosteric effectors. These do not need to be included in the reaction's pattern file.
 - ``allosteric``: whether or not the reaction is allosteric. In general you set this to ``1`` if there are negative/positive effectors, otherwise it's set to ``0``.
 - ``subunits``: how many subunits is the enzyme composed by. Mostly relevant for allosteric reactions.



.. _pattern_file_specification:

Pattern file specification
-----------------------------------

To define the mechanism for a given reaction/enzyme, GRASP needs a pattern file to be specified. These files should be inside the ``patterns`` folder, where you can find some examples as well.

Here we will show how to specify a pattern file.

Let's consider a uni-uni mechanism as the one below:
::

    E_c + m_3pg_c <-> E_c&m_3pg_c
    E_c&m_3pg_c <-> E_c&m_2pg_c
    E_c&m_2pg_c <-> E_c + m_2pg_c


where each line represents one elementary reaction, ``E`` is the enzyme, ``3pg`` is the substrate, ``2pg`` is the product, ``E_c&m_3pg_c`` is the complex of the enzyme ``E`` and the substrate ``3pg``, and ``E_c&m_2pg_c`` is the complex of the enzyme ``E`` and the substrate ``2pg``. ``_c`` denotes the compartment where the metabolites and enzyme are located, and it is optional.

The pattern file that needs to be specified for this reaction mechanism would look like:
::

    1 2 k01.*A
    2 1 k02
    2 3 k03
    3 2 k04
    3 1 k05
    1 3 k06.*P

where ``k01, ..., k06`` are the rate constants for the elementary reactions, ``A`` is the substrate, and ``P`` the product. |br|
The two numbers at the beginning of each row represent the different enzyme states: the free enzyme (1), the enzyme bound to ``A`` (2), and the enzyme bound to ``B`` (3). |br|
The first of the two numbers is the starting state and the second is the end state, e.g. in the first line ``1 2 k01.*A``, ``1`` is the starting state, the free enzyme, and ``2`` is the end state, the enzyme bound to substrate ``A``. |br|
Since all elementary reactions must be reversible, two lines are specified for each elementary reaction, one for the forward direction, e.g. ``1 2 k01.*A`` where ``A`` binds to the free enzyme, and another for the reverse direction, e.g. ``2 1 k02`` where ``A`` is released.

The convention for metabolite names is:

  - A, B, C, D refer to substrates;
  - P, Q, R, S refer to products;
  - I refers to inhibitors; if there is more than one, try I1, I2, I3 and make sure it worked by looking into the reactions functions;
  - Z refers to the activators; if there is more than one, try Z1, Z2, Z3 and make sure it worked by looking into the reactions functions.

You can also use the package `set_up_grasp_models <https://github.com/biosustain/set_up_grasp_models>`_ to generate the pattern files from a file with the elementary reactions.



Example
-----------------------------------


To build a model ensemble you can use the script ``build_model.m`` in the examples folder, which looks similar to the code below:


.. code-block:: matlab

    % maximum number of models sampled, no matter what
    maxNumberOfSamples = 10000;

    % threshold of the jacobian's eigenvalues
    eigThreshold = 10^-5;

    modelID = 'toy_model';
    inputFile = fullfile('..', 'io', 'input', modelID);
    outputFile = fullfile('..', 'io','output', [modelID, '.mat']);

    ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);


To build the model ensemble we use the ``buildEnsemble`` function, which takes 4 arguments:

 - ``inputFile``: path to the input excel file;
 - ``outpuFile``: path to the output file (to be created);
 - ``maxNumberOfSamples``: the maximum number of models that will be sampled no matter how many models are valid. The goal is to get the number of valid models specified in the excel input file. However, if there are no valid models (very unlikely) the program will run forever. ``maxNumberOfSamples`` will prevent that. A model is invalid if the real part of the eigenvalues of its jacobian is higher than ``eigThreshold`;
 - ``eigThreshold``: threshold for the real part of the model's jacobian eigenvalues. Models with a eigenvalue real part higher than ``eigThreshold`` are discarded.

In general it is recommended to use the ``io/input`` and ``io/output`` folders to store your input/output files, but you can use any other folders as long as the paths are specified correctly.


.. |br| raw:: html

      <br>