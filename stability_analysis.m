clear				            																
addpath('./patternFxns','./ensembleFxns');


% threshold of the jacobian's eigenvalues
eigThreshold = 10^-5;
modelID = 'HMP2360_r0_t3_no_promiscuous2';
outputFolder = 'output_test/';


load(strcat(outputFolder, modelID, '.mat'));

stabilityRes = ensembleStabilityTest(ensemble, eigThreshold);

save(strcat(outputFolder, 'stability_', modelID, '.mat'), 'stabilityRes');