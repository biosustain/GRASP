classdef sampleGibbsReactionEnergiesTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testSampleGibbsReactionEnergies1(testCase)
            
            seed = 1;
            rng(seed)
           
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_measuredMets'));
            ensemble = ensemble.ensemble;
            models(1).poolFactor = [];
            strucIdx = 1;
 
            [ensemble, models] = sampleGibbsReactionEnergies(ensemble, models, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleGibbsReactionEnergies1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
                   
            testCase.verifyEqual(trueResModels, models);
            testCase.verifyEqual(trueResEnsemble, ensemble);
        end
       
    end
end


