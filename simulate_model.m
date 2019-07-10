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


% If you have reference metabolite concentrations provide them below,
% otherwise comment the code and remove metsRefConc from the argument
% list.
inputFolder = 'input_test';
metsRefConc = readtable(fullfile(inputFolder, 'ref_mets_v2_3.csv'));
metsRefConc = table2array(metsRefConc(:,2));
metsRefConc = metsRefConc(ensemble.metsActive);


% Change initial conditions here if you want
enzymesIC(2) = 1.5;
metsIC(5) = 2;

% Specifiy the time of simulation (probably in hours)
finalTime = 1;

simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsRefConc);

save(fullfile(outputFolder, ['simulation_', modelID, '.mat']), 'simulationRes')
write(cell2table(ensemble.mets(ensemble.metsActive)), fullfile(outputFolder, [modelID, '_metsActive.dat']));
write(cell2table(ensemble.rxns(ensemble.activeRxns)), fullfile(outputFolder, [modelID, '_rxnsActive.dat']));



