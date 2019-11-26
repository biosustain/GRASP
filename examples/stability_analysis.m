% Calculate the jacobian eigenvalues and check if the real part is larger
% than the given threshold

clear				            																
addpath(fullfile('..', 'matlab_code', 'analysisFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'), ...
        fullfile('..', 'matlab_code', 'patternFxns'));

% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;
modelID = 'toy_model';
outputFolder = fullfile('..', 'io','output');


load(fullfile(outputFolder, [modelID, '.mat']));

stabilityRes = ensembleStabilityTest(ensemble, eigThreshold);

save(fullfile(outputFolder, ['stability_', modelID, '.mat']), 'stabilityRes');