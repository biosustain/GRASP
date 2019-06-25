% Simulate a model ensemble

clear
addpath('patternFxns','ensembleFxns');


modelID = 'HMP2360_r0_t3_new';
outputFolder = 'output_test';

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

simulationRes = simulateEnsemble(ensemble, enzymesIC, metsIC);

save(fullfile(outputFolder, ['simulation_', modelID, '.mat']), 'simulationRes')


