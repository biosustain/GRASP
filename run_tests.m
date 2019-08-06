%
% Test cases for GRASP
%
%% Case 1: Nick's model Glycolysis_Grasp_isoenzymes
clear
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/glycolysis_isoenzymes');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp([newline, newline, '--- glycolysis_isoenzymes ran fine ---',newline, newline]);

%% Case 2: Pedro's MEP_example
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

disp([newline, newline, '--- MEP_example ran fine ---',newline, newline]);


%% Case 3: Pedro's MEP_example_inhibitors, to test inhibitors
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

disp([newline, newline, '--- MEP_example_inhibitors ran fine ---', newline, newline]);


%% Case 4: HMP2360_r0_t0_nick, all models should be set up without issues

clearvars
rng('shuffle');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('patternFxns','ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/HMP2360_r0_t0_nick');           % Here the test case HMP pathway model is chosen

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);

disp([newline, newline, '--- HMP2360_r0_t0_nick ran fine ---', newline, newline]);


%% Case 5: HMP2360_r0_t3_new, all models should be stable

clearvars
rng('shuffle');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('patternFxns','ensembleFxns');

% only valid/stable models are kept, and it will keep sampling until the
%  "Number of particles" defined in the excel is reached, however, it is a
%  good idea to define the maximum number of models to be sampled in total, 
%  otherwise if no stable models are found it will go on sampling forever.
maxNumberOfSamples = 10000;   

% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;

modelID = 'HMP2360_r0_t3_new';
inputFile = fullfile('input_test', modelID);
outputFile = fullfile('output_test', [modelID, '.mat']);

ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);


disp([newline, newline, '--- HMP2360_r0_t3_new ran fine ---', newline, newline]);


%% Case 6: HMP2360_r0_t3_no_promiscuous2, includes random mechanisms, about 10-20% of models were usually unstable before the Gibbs energy fix, now less than 5% should be unstable.

clearvars
rng('shuffle');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('patternFxns','ensembleFxns');

% only valid/stable models are kept, and it will keep sampling until the
%  "Number of particles" defined in the excel is reached, however, it is a
%  good idea to define the maximum number of models to be sampled in total, 
%  otherwise if no stable models are found it will go on sampling forever.
maxNumberOfSamples = 10000;   

% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;

modelID = 'HMP2360_r0_t3_no_promiscuous2';
inputFile = fullfile('input_test', modelID);
outputFile = fullfile('output_test', [modelID, '.mat']);

ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);


disp([newline, newline,'--- HMP2360_r0_t3_no_promiscuous2 ran fine ---', newline, newline]);


%% Case 7: HMP1489_r1_t0_test_gibbs, should crash with an error on TMFA

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

disp([newline, newline, '--- HMP1489_r1_t0_test_gibbs ran fine ---', newline, newline]);