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
            
            loadEnsembleStructure(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2_sbml'));
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2_sbml.mat'));
            ensemble = ensemble.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            modelI = 1;
            exportSBML(ensemble, modelI, outputFolder);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_no_promiscuous2_1.xml'));
            SBMLres = fileread(filepath);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_toy_model1_no_promiscuous2_1.xml'));
            trueRes = fileread(filepath);
            
            % Need to skip Simbiology version in the xml file, otherwise
            % the test will fail for different versions.
            indTrueRes = strfind(trueRes, '<model id=');
            indRes = strfind(SBMLres, '<model id=');
            testCase.verifyEqual(SBMLres(indRes:end), trueRes(indTrueRes:end));
        end
       
        
       function testExportSBMLPromiscuous(testCase)
            
            seed = 1;
            rng(seed)
            
            loadEnsembleStructure(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2_sbml'));
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2_sbml.mat'));
            ensemble = ensemble.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            modelI = 1;
            exportSBML(ensemble, modelI, outputFolder);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric2_1.xml'));
            SBMLres = fileread(filepath);
            
            filepath = (fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_toy_model1_allosteric2_1.xml'));
            trueRes = fileread(filepath);
            
            % Need to skip Simbiology version in the xml file, otherwise
            % the test will fail for different versions.       
            indTrueRes = strfind(trueRes, '<model id=');
            indRes = strfind(SBMLres, '<model id=');
            testCase.verifyEqual(SBMLres(indRes:end), trueRes(indTrueRes:end));
        end
        
    end
end

