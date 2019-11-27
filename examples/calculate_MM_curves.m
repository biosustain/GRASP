% Calculate the MM curves for each reaction-metabolite pair.

clearvars
      				            																% check first that no other process is running
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));
    

modelID = 'toy_model';
outputFolder = fullfile('..', 'io','output');

numModels = 5;

load(fullfile(outputFolder, [modelID, '.mat']));

% Optional arguments
structIdx = 1;  % model structure ID
saturatingConc = 10^2*10^6;  %saturating concentration
substrateRange = logspace(-15,6);  % range of substrate concentrations for which fluxes will be calculated
rxnList = 1:3;  % reactions for which the MM curves will be calculated

% Add reactions files to path
reactionsFolder = fullfile('..', 'reactions', [modelID, '_', num2str(structIdx)]);
addpath(reactionsFolder);

calculateMMCurves(outputFolder, ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
