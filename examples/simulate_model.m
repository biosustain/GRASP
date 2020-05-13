% Simulate a model ensemble

clear
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

modelID = 'toy_model';
outputFolder = fullfile('..', 'io','output');

% How many models in the ensemble you want to simulate
numModels = 1000;

% Define how many seconds until ODE solver is interrupted. The idea is to
% skip models that take ages to simulate.
interruptTime = 40;

load(fullfile(outputFolder, [modelID, '.mat']))

% Define whether the initial condition for metabolites is a relative or an
% absolute concentration by setting metsAbsOrRel to either 'rel' or 'abs',
% respectively
metsAbsOrRel = 'rel';

% Change initial conditions here if you want, format: {rxn/met ID, initial value}
enzymesIC = {{'r1', 0.2}, {'r10', 2}};       % Always relative concentrations for enzymes
metsIC = {{'m5', 0.6}, {'m14', 1.5}};        % Absolute concentrations must be given in mol/L

% Specifiy the time of simulation (probably in hours)
finalTime = 1;

simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime, numModels);

save(fullfile(outputFolder, ['simulation_', modelID, '.mat']), 'simulationRes')
write(cell2table(ensemble.mets(ensemble.metsActive)), fullfile(outputFolder, [modelID, '_metsActive.dat']));
write(cell2table(ensemble.rxns(ensemble.activeRxns)), fullfile(outputFolder, [modelID, '_rxnsActive.dat']));



