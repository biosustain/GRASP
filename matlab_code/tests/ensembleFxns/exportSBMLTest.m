classdef exportSBMLTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testExportSBMLNoPromiscuous(testCase)
            
            seed = 1;
            rng(seed)
            
            loadEnsembleStructure(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2'))
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2.mat'));
            ensemble = ensemble.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            modelI = 1;
            exportSBML(ensemble, modelI, outputFolder);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2_1.xml'));
            SBMLres = fileread(filepath);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_toy_model1_no_promiscuous2_1.xml'));
            trueRes = fileread(filepath);
                   
            testCase.verifyEqual(trueRes, SBMLres);
        end
       
        
       function testExportSBMLPromiscuous(testCase)
            
            seed = 1;
            rng(seed)
            
            loadEnsembleStructure(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2'))
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2.mat'));
            ensemble = ensemble.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            modelI = 1;
            exportSBML(ensemble, modelI, outputFolder);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2_1.xml'));
            SBMLres = fileread(filepath);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_toy_model1_allosteric2_1.xml'));
            trueRes = fileread(filepath);
                   
            testCase.verifyEqual(trueRes, SBMLres);
        end
        
    end
end

