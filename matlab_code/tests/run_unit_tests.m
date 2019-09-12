clearvars

import matlab.unittest.TestSuite;
suitePatterFxns = TestSuite.fromFolder('patternFxns');

						   				            																% check first that no other process is running
addpath(fullfile('..',  'patternFxns'), ...
        fullfile('..',  'ensembleFxns'), ...
        'ensembleFxns');
    
    
currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
tempReactionsFolder = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions');
mkdir(tempReactionsFolder);


result = run(suitePatterFxns)

% ensembleFxns tests - temporary
run(loadEnsembleStructureTest)
run(initializeEnsembleTest)
run(initialSamplerTest)
run(sampleEnzymeAbundancesTest)



if exist(tempReactionsFolder, 'dir')
    rmdir(tempReactionsFolder, 's');
end