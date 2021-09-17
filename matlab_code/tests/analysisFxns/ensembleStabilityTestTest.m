classdef ensembleStabilityTestTest < matlab.unittest.TestCase
   
    properties
        currentPath
        relTol = 1e-10;
        absTol = 1e-10;
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1');
            
            if exist(reactionsFolder, 'dir')
                rmdir(reactionsFolder, 's');
            end
        end
    end
    
    
    methods (Test)
        function testEnsembleStabilityTest1(testCase)
                        
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            ensemble.freeVars{end+1} = 'r_r13';
            
            eigThreshold = -0.5;
            stabilityRes = ensembleStabilityTest(ensemble,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResEnsembleStabilityTest1'));
            trueRes = trueRes.stabilityRes;
                   
            testCase.verifyThat(stabilityRes, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
    end
end

