% Simulate a model ensemble

clear
addpath(fullfile('..', 'matlab_code', 'patternFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'));


modelID = 'HMP2360_r0_t3_new';
outputFolder = 'output_test';

% Define how many seconds until ODE solver is interrupted. The idea is to
% skip models that take ages to simulate.
interruptTime = 40;

load(fullfile(outputFolder, [modelID, '.mat']))


% Get default initial conditions, all ones
freeVars = numel(ensemble.freeVars);
xopt = ones(freeVars,1);
ix_mets = 1:numel(ensemble.metsActive);
ix_enz = ix_mets(end)+1:freeVars;
metsIC = xopt(ix_mets);
enzymesIC = xopt(ix_enz);


% Change initial conditions here if you want
enzymesIC(2) = 1.5;
metsIC(5) = 2;

% Specifiy the time of simulation (probably in hours)
finalTime = 1;

simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, interruptTime);

save(fullfile(outputFolder, ['simulation_', modelID, '.mat']), 'simulationRes')
write(cell2table(ensemble.mets(ensemble.metsActive)), fullfile(outputFolder, [modelID, '_metsActive.dat']));
write(cell2table(ensemble.rxns(ensemble.activeRxns)), fullfile(outputFolder, [modelID, '_rxnsActive.dat']));


