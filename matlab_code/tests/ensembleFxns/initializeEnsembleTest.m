classdef initializeEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
    
 
    methods (Test)
        function testInitializeEnsemble(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleWithRandom(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_random.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
    end
end

