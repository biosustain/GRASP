classdef sampleFluxesAndGibbsFreeEnergiesTest < matlab.unittest.TestCase

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
        function testSampleFluxesAndGibbsFreeEnergiesNormal(testCase)
            
            seed = 1;
            rng(seed)

            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs.mat'));
            ensemble = ensemble.ensemble;
            maxNumberOfSamples = 100;
            ensemble.fluxPrior = 'normal';
            ensemble.thermoPrior = 'normal';
            
            ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleFluxesAndGibbsFreeEnergiesNormal.mat'));
            trueResEnsemble = trueRes.ensemble;
           
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));    
        end
        
        function testSampleFluxesAndGibbsFreeEnergiesUniform(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs'));
            ensemble = ensemble.ensemble;
            maxNumberOfSamples = 100;
            ensemble.fluxPrior = 'uniform';
            ensemble.thermoPrior = 'uniform';
            
            ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleFluxesAndGibbsFreeEnergiesUniform.mat'));
            trueResEnsemble = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end


