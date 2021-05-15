classdef controlAndResponseAnalysisTest < matlab.unittest.TestCase
   
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
        function testControlAndResponseAnalysis1(testCase)
            
            
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            ensemble.freeVars{end+1} = 'r_r13';
            
            saveResMatrices = 0;
            mcaResults = controlAndResponseAnalysis(ensemble,saveResMatrices);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResControlAndResponseAnalysis1'));
            trueRes = trueRes.mcaResults;
                   
            testCase.verifyThat(mcaResults, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testControlAndResponseAnalysisNumbers(testCase)
                        
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_numbers');
            loadEnsembleStructure(xlsxFile);

            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_numbers.mat'));
            ensemble = ensemble.ensemble;
            
            saveResMatrices = 0;
            mcaResults = controlAndResponseAnalysis(ensemble,saveResMatrices);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResControlAndResponseAnalysisNumbers'));
            trueRes = trueRes.mcaResults;

            testCase.verifyThat(mcaResults, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end

