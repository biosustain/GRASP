classdef checkStabilityTest < matlab.unittest.TestCase
   
    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2');
            
            if exist(reactionsFolder, 'dir')
                rmdir(reactionsFolder, 's');
            end
        end
    end
    
    
    methods (Test)
        function testCheckStability1(testCase)
            
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            models = load(fullfile(testCase.currentPath{1}, 'testFiles', 'models_toy_model1_random2_checkStab.mat'));
            models = models.models;
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1_random2_checkStab.mat'));
            ensemble = ensemble.ensemble;
            ensemble.freeVars{end+1} = 'r_r13';
            [~,ensemble.metsLI] = rref(ensemble.Sred');
            
            strucIdx = 1;
            eigThreshold = 10^-5;
            isModelStable = checkStability(ensemble,models,strucIdx,eigThreshold);
            
            trueRes = true;
                   
            testCase.verifyEqual(isModelStable, trueRes)
        end
    end
end

