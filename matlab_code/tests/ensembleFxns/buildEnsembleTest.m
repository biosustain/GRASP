classdef buildEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_allosteric2_1')};
            
            for i = 1:size(reactionsFolderList, 1)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
    end
    
    
    methods (Test)
        function testBuildEnsembleRandom(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_random2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleRandom'));
            trueRes = trueRes.ensemble;
                   
            testCase.verifyEqual(trueRes, ensemble);
        end
       
        function testBuildEnsembleAllosteric(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_allosteric2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleAllosteric'));
            trueRes = trueRes.ensemble;
                   
            testCase.verifyEqual(trueRes, ensemble);
        end
       
    end
end




