% Sample a kinetic model ensemble (reference point)

clearvars
rng('shuffle');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath(fullfile('..', 'matlab_code', 'patternFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'));

% only valid/stable models are kept, and it will keep sampling until the
%  "Number of particles" defined in the excel is reached, however, it is a
%  good idea to define the maximum number of models to be sampled in total, 
%  otherwise if no stable models are found it will go on sampling forever.
maxNumberOfSamples = 10000;   

% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;

modelID = 'toy_model';
inputFile = fullfile('..', 'io', 'input', modelID);
outputFile = fullfile('..', 'io','output', [modelID, '.mat']);

ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold);
