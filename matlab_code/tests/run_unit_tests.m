clearvars

import matlab.unittest.TestSuite;
suitePatterFxns = TestSuite.fromFolder('patternFxns');

						   				            																% check first that no other process is running
addpath(fullfile('..',  'patternFxns'), ...
        fullfile('..',  'ensembleFxns'));
    
    
currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
tempReactionsFolder = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions');
mkdir(tempReactionsFolder);
            
%testCase = loadEnsembleStructureTest;
%testCase = initializeEnsembleTest;
%testCase = initialSamplerTest;
%testCase = readInputTest;
%res = run(testCase)

result = run(suitePatterFxns)

if exist(tempReactionsFolder, 'dir')
    rmdir(tempReactionsFolder, 's');
end