Tutorials
=================


Building a model ensemble
--------------------------------

To build a model ensemble you can use the script `build_model.m` in the examples folder, which looks similar to the code below:


.. code-block:: matlab

    % maximum number of models sampled, no matter what
    maxNumberOfSamples = 10000;   

    % threshold of the jacobian's eigenvalues
    eigThreshold = 10^-5;

    modelID = 'toy_model';
    inputFile = fullfile('..', 'io', 'input', modelID);
    outputFile = fullfile('..', 'io','output', [modelID, '.mat']);

    ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);


To build the model ensemble we use the `buildEnsemble` function, which takes 4 arguments:

 - `inputFile`: path to the input file;
 - `outpuFile`: path to the output file;
 - `maxNumberOfSamples`: the maximum number of models that will be sampled no matter how many models are valid. The goal is to get the number of valid models specified in the excel input file, however, if there are no valid models (very unlikely) the program will run forever, `maxNumberOfSamples` will prevent that;
 - `eigThreshold`: threshold for the real part of the model's jacobian eigenvalues. Models with a eigenvalue real part higher than `eigThreshold` are discarded.

In general it is recommended to use the `io/input` and `io/output` folders to store your input/output files, but you can use any other folders as long as you specify the path correctly.



Setting up the input excel file
--------------------------------

To build a model ensemble you need to specify the model first, this is done by setting up an input excel file where the following is specified:

  - reactions stoichiometry;
  - reaction fluxes (in micromol/gCDW/h);
  - metabolite concentrations (in mol/L) for thermodynamics to calculate the reactions Gibbs energies;
  - Gibbs standard energies;
  - enzyme mechanisms;
  - enzyme regulation.

In the folder `io/input` you can find an example.
However, since this is a tedious and error prone process, we've created a python package that automates most of the process and validates the final file, see `here <https://github.com/biosustain/set_up_grasp_models>`_. 

For a more thorough description of the input file, please see how_to_build_a_model_.


Setting up the mechanism files
--------------------------------

To define the mechanism for a given reaction/enzyme, GRASP needs a pattern file to be specified.
For a uni-uni mechanism as the one below:
::

    E_c + m_3pg_c <-> E_c&m_3pg_c
    E_c&m_3pg_c <-> E_c&m_2pg_c
    E_c&m_2pg_c <-> E_c + m_2pg_c


where each line represents one elementary reaction, E is the enzyme, 3pg is the substrate, and 2pg is the product. 

The pattern file that needs to be specified would look like:
::
    1 2 k01.*A
    2 1 k02
    2 3 k03
    3 2 k04
    3 1 k05
    1 3 k06.*P

where k01, ..., k06 are the rate constants for the elementary reactions, A is the substrate, and P the product. The two numbers at the beginning of each row represent the different enzyme states: the free enzyme (1), the enzyme bound to A (2), and the enzyme bound to B (3). The first of the two numbers is the starting state and the second is the end state, e.g. in the first line `1 2 k01.*A`, `1` is the starting state, the free enzyme, and `2` is the end state, the enzyme bound to substrate A.
Since all elementary reactions must be reversible, two lines are specified for each elementary reaction, one for the forward direction, e.g. `1 2 k01.*A` where A binds to the free enzyme, and another for the reverse direction, e.g. `2 1 k02` where A is released.
The convention for metabolite names is:
When specifying pattern files, the convention is:

  - A, B, C, D refer to substrates;
  - P, Q, R, S refer to products;
  - I refers to inhibitors, if there is more than one, try I1, I2, I3 and make sure it worked by looking into the reactions functions;
  -  Z refers to the activators, if there is more than one, try Z1, Z2, Z3 and make sure it worked by looking into the reactions functions.

You can also use the package `set_up_grasp_models <https://github.com/biosustain/set_up_grasp_models>`_ to generate the pattern files from a file with the elementary reactions.


Simulating the model ensemble
--------------------------------



Analyzing the model ensemble
--------------------------------


Metabolic Control Analysis (MCA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stability analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Visualization (with python)
--------------------------------