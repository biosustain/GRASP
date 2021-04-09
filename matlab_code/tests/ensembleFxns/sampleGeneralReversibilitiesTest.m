classdef sampleGeneralReversibilitiesTest < matlab.unittest.TestCase

    properties
        currentPath
        relTol = 1e-4;
        absTol = 1e-4;
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testSampleGeneralReversibilities(testCase)
            
            seed = 1;
            rng(seed)
           
            load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampleGeneralReversibilities_toy_model1'));

            [ensemble, models, isModelValid] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleGeneralReversibilities1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
            
            testCase.verifyEqual(isModelValid, true);
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueResModels, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testSampleGeneralReversibilitiesRandom(testCase)
            
            seed = 1;
            rng(seed)
           
            load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampleGeneralReversibilities_toy_model_random2'));
            

            [ensemble, models, isModelValid] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleGeneralReversibilitiesRandom'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
            
            testCase.verifyEqual(isModelValid, true);
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueResModels, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end


