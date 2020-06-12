% Simulate a model ensemble

% Clear all variables and add functions to path
clear
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

% Define the model ID, which should basically be the filename of the model
%  ensemble
modelID = 'toy_model';
outputFolder = fullfile('..', 'io','output');

% How many models in the ensemble you want to simulate
numModels = 5;

% How many cores you want to use to run the simulation
numCores = 2;

% Define how many seconds until ODE solver is interrupted. The idea is to
% skip models that take too long to simulate.
interruptTime = 40;

% Load model ensemble
load(fullfile(outputFolder, [modelID, '.mat']))

% Define whether the initial condition for metabolites is a relative or an
% absolute concentration by setting metsAbsOrRel to either 'rel' or 'abs',
% respectively
metsAbsOrRel = 'rel';

% Change initial conditions here if you want, format: {rxn/met ID, initial value}
enzymesIC = {{'r1', 1}, {'r10', 1}};       % Always relative concentrations for enzymes
metsIC = {{'m5', 0.6}, {'m11', 1.2}};      % When setting metsAbsOrRel = 'abs', absolute concentrations must be given in mol/L

% Specifiy the time of simulation according to the flux units you've 
%  provided in the input excel file (probably in hours)
finalTime = 1;


% Fluxes will be in absolute units, while metabolite concentrations will be
%  scaled.
simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime, numModels, numCores);

% Save the results in a .mat file that can be imported in python for 
%  further analysis. This might not work for very large files.
save(fullfile(outputFolder, ['simulation_', modelID, '.mat']), 'simulationRes')
write(cell2table(ensemble.mets(ensemble.metsActive)), fullfile(outputFolder, [modelID, '_metsActive.dat']));
write(cell2table(ensemble.rxns(ensemble.activeRxns)), fullfile(outputFolder, [modelID, '_rxnsActive.dat']));



