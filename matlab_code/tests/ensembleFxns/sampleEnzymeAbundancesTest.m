classdef sampleEnzymeAbundancesTest < matlab.unittest.TestCase

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
        function testSampleEnzymeAbundances1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1'));
            ensemble = ensemble.ensemble;
            models(1).poolFactor = [];
            strucIdx = 1;
            
            [models] = sampleEnzymeAbundances(ensemble,models,strucIdx);
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleEnzymeAbundances1'));
            trueRes = trueRes.models;
                   
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testSampleEnzymeAbundancesRandom(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_random'));
            ensemble = ensemble.ensemble;
            models(1).poolFactor = [];
            strucIdx = 1;
            
            [models] = sampleEnzymeAbundances(ensemble,models,strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleEnzymeAbundancesRandom'));
            trueRes = trueRes.models;
                   
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end

