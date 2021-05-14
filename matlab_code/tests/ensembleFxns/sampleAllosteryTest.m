classdef sampleAllosteryTest < matlab.unittest.TestCase

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
        function testSampleAllostery1(testCase)
            
            seed = 1;
            rng(seed)

            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampleAllostery_toy_model1_allosteric2.mat'));
            ensemble = ensemble.ensemble;
            models = ensemble.populations.models;
            
            for i=1:numel(models)
                models(i).rxnParams = rmfield(models(i).rxnParams, 'L');
                models(i).rxnParams = rmfield(models(i).rxnParams, 'allostericFactors');
                models(i).rxnParams = rmfield(models(i).rxnParams, 'KposEff');
                models(i).rxnParams = rmfield(models(i).rxnParams, 'KnegEff');
            end
                                    
            strucIdx = 1;
            [ensemble, models] = sampleAllostery(ensemble, models, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleAllosteryTest1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
                   
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueResModels, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
       
    end
end

