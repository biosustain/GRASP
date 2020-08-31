classdef controlAnalysisTest < matlab.unittest.TestCase
   
    properties
        currentPath
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
        function testControlAnalysis1(testCase)
                        
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            ensemble.freeVars{end+1} = 'r_r13';
            
            saveResMatrices = 0;
            mcaResults = controlAnalysis(ensemble,saveResMatrices);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResControlAnalysis1'));
            trueRes = trueRes.mcaResults;

            testCase.verifyThat(mcaResults, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.AbsoluteTolerance(1e-10)));
        end
    end
end

