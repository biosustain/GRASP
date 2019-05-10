% Example 1. Sample a kinetic model of the MEP pathway (reference point)
%--------------------------------------------------------------------------
% Executes GRASP workflow
%
% Inputs:       (-)
%
% Outputs:      (-)
%--------------------- Pedro Saa 2017 -------------------------------------
clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/Glycolysis_Grasp_isoenzymes');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp('GRASP_input_yeast ran fine');

clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/MEP_example');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp('MEP_example ran fine');


clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/MEP_example_inhibitors');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp('MEP_example_inhibitors ran fine');




clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/HMP2360_r0_t0_nick');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp('HMP2360_r0_t0_nick ran fine');


% This last tests throws an error on purpose, the error message should be:
% "The TMFA problem is infeasible. Verify that the standard Gibbs free energy and
% metabolite concentration values are valid/correct. Reactions
% r_AANAT, r_AANAT_tryptm with standard Gibbs energies'[100,200][10,20]' seem to be the problem."

clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/HMP1489_r1_t0_test_gibbs');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp('HMP1489_r1_t0_test_gibbs ran fine');