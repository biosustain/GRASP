classdef loadEnsembleStructureTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_allosteric_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_freeExchange_1')};
            
            for i = 1:size(reactionsFolderList, 1)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
    end
       
    
    methods (Test)
        function testLoadEnsembleStructure1(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructure1WithAllostery(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_allosteric.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructure1WithRandom(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_random.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructureWithFreeExchange(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_freeExchange');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_freeExchange.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
    end
end

