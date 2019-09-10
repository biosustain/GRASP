% Run MCA analysis on your model ensemble

clear      				   
addpath(fullfile('..', 'matlab_code', 'patternFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'));

saveMCAMatrices = 1;    % whether or not to save the MCA results for all models and not just mean values
modelID = 'HMP2360_r0_t3_new';
outputFolder = 'output_test';

load(fullfile(outputFolder, [modelID, '.mat']))


% Run MCA analysis

mcaResults = controlAnalysis(ensemble, saveMCAMatrices);

% If you have promiscuous enzymes or isoenzymes you might want to run
% controlAndResponseAnalysis instead
% mcaResults = controlAndResponseAnalysis(ensemble, saveMCAMatrices);


% Save MCA results
save(fullfile(outputFolder, ['MCA_', modelID, '.mat']), 'mcaResults');
write(cell2table(ensemble.rxns(ensemble.activeRxns)), fullfile(outputFolder, [modelID, '_rxnsActive.dat']));
write(cell2table(ensemble.mets(ensemble.metsActive)), fullfile(outputFolder, [modelID, '_metsActive.dat']));
write(cell2table(mcaResults.enzNames), fullfile(outputFolder, [modelID, '_enzNames.dat']));


% Plot MCA results

% Optional, Define ranges for displaying the MCA results: {1st category, range; 2nd category, range}  
% For example, categories = {'Glycolysis',[1,20]; 'Pentose Phosphate Pathway',[21,30];'Others', [31,37] ;'Ethanol metabolism', [38,49];'TCA cycle', [50,76] ;'ATP ADP NADH', [77,83]; 'Exchange reactions', [84,88]};
categories = {'all', [1, 5]}; % currently not working when empty

plotControlAnalysis(mcaResults, ensemble, categories);

% If doing controlAndResponseAnalysis choose plotControlAndResponseAnalysis
% instead
%plotControlAndResponseAnalysis(mcaResults, ensemble);


