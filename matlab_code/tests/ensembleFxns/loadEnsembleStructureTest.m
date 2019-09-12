classdef loadEnsembleStructureTest < matlab.unittest.TestCase
   
    methods (Test)
        function testLoadEnsembleStructure1(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructure1WithAllostery(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1_allosteric');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_allosteric.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructure1WithRandom(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1_random');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_random.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
        
        function testLoadEnsembleStructureWithFreeExchange(testCase)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1_freeExchange');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_freeExchange.mat'));
            
            testCase.verifyEqual(trueRes.ensemble,ensemble);            
        end
    end
end

