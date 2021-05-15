classdef plotControlAndResponseAnalysisTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testPlotControlAndResponseAnalysisTest1(testCase)
                        
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            
            mcaResults = load(fullfile(testCase.currentPath{1}, 'testFiles', 'MCA_response_toy_model1_random2.mat'));
            mcaResults = mcaResults.mcaResults;            
            
            categories = {};

            plotControlAndResponseAnalysis(mcaResults, ensemble);

        end
        
        function testPlotControlAndResponseAnalysisTestCategories1(testCase)
                        
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            
            mcaResults = load(fullfile(testCase.currentPath{1}, 'testFiles', 'MCA_response_toy_model1_random2.mat'));
            mcaResults = mcaResults.mcaResults;
                        
            categories = {'all', [1, 5]};
            enzymeCategories = {'all', [1, 4]};

            plotControlAndResponseAnalysis(mcaResults, ensemble, categories, enzymeCategories);

        end
        
        function testPlotControlAndResponseAnalysisTestNumbers(testCase)
                        
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_numbers.mat'));
            ensemble = ensemble.ensemble;
            
            mcaResults = load(fullfile(testCase.currentPath{1}, 'testFiles', 'mca_response_toy_model1_numbers.mat'));
            mcaResults = mcaResults.mcaResults;            
            
            categories = {};
            
            plotControlAnalysis(mcaResults, ensemble, categories);

        end
    end
end