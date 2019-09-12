% Calculate the jacobian eigenvalues and check if the real part is larger
% than the given threshold

clear				            																
addpath(fullfile('..', 'matlab_code', 'patternFxns'), ...
        fullfile('..', 'matlab_code', 'ensembleFxns'));

% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;
modelID = 'HMP2360_r0_t3_no_promiscuous2';
outputFolder = 'output_test';


load(fullfile(outputFolder, [modelID, '.mat']));

stabilityRes = ensembleStabilityTest(ensemble, eigThreshold);

save(fullfile(outputFolder, ['stability_', modelID, '.mat']), 'stabilityRes');