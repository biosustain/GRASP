classdef sampleModifierElemFluxesTest < matlab.unittest.TestCase

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
        function testSampleModifierElemFluxes1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_withGibbsRever_toy_model1_allosteric2'));
            ensemble = ensemble.ensemble;
            
            models = load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampledEnzymeAbundancesGibbsRever_toy_model1_allosteric2'));
            models = models.models;
            
            strucIdx = 1;
            models = sampleModifierElemFluxes(ensemble, models, strucIdx);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleModifierElemFluxes1'));
            trueRes= trueRes.models;
                   
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testSampleModifierElemFluxesRandom(testCase)
            
            seed = 1;
            rng(seed)
                       
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_withGibbsRever_toy_model1_random2'));
            ensemble = ensemble.ensemble;
            
            models = load(fullfile(testCase.currentPath{1}, 'testFiles', 'sampledEnzymeAbundancesGibbsRever_toy_model1_random2'));
            models = models.models;
            
            strucIdx = 1;
            models = sampleModifierElemFluxes(ensemble, models, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleModifierElemFluxesRandom'));
            trueRes= trueRes.models;
                   
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  

        end
       
    end
end
