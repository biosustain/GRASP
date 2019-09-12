classdef initializeEnsembleTest < matlab.unittest.TestCase

    methods (Test)
        function testInitializeEnsemble(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            filepath = fullfile(currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleWithRandom(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            filepath = fullfile(currentPath{1}, 'testFiles', 'ensemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_random.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
    end
end

