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
        function testSampleGibbsReactionEnergiesGurobi(testCase)
            
            seed = 1;
            rng(seed)
           
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_measuredMets'));
            ensemble = ensemble.ensemble;
            models(1).poolFactor = [];
            strucIdx = 1;
            
            solver  = 'gurobi';
            ensemble.LPSolver = solver;            
 
            [ensemble, models] = sampleGibbsReactionEnergies(ensemble, models, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleGibbsReactionEnergies1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
            trueResEnsemble.LPSolver = solver;
           
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueResModels, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)));
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)));
        end
       
        function testSampleGibbsReactionEnergiesLinprog(testCase)
            
            seed = 1;
            rng(seed)
           
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_measuredMets'));
            ensemble = ensemble.ensemble;
            models(1).poolFactor = [];
            strucIdx = 1;
            
            solver  = 'linprog';
            ensemble.LPSolver = solver;            
 
            [ensemble, models] = sampleGibbsReactionEnergies(ensemble, models, strucIdx);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSampleGibbsReactionEnergies1'));
            trueResModels = trueRes.models;
            trueResEnsemble = trueRes.ensemble;
            trueResEnsemble.LPSolver = solver;
           
            testCase.verifyThat(models, matlab.unittest.constraints.IsEqualTo(trueResModels, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)));
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueResEnsemble, ...
                 'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)));
        end
    end
end


