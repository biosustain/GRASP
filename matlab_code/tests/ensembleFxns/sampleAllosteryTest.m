classdef sampleAllosteryTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testSampleAllostery1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_allosteric2'));
            ensemble = ensemble.ensemble;
            
            models = load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampledEnzymeAbundances_toy_model1_allosteric2'));
            models = models.models;
            
            strucIdx = 1;
            [ensemble, models] = sampleAllostery(ensemble, models, strucIdx);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleAllosteryTest1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
                   
            testCase.verifyEqual(trueResModels, models);
            testCase.verifyEqual(trueResEnsemble, ensemble);
        end
       
    end
end

