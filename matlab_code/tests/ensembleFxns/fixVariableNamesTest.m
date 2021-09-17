classdef fixVariableNamesTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testFixVariableNames1(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            [xRxns,rxnsList] = xlsread(xlsxFile,'rxns');                         % load rxn info            
           
            rxnsList = fixVariableNames(rxnsList, 'r');

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResFixVariableNames1'));
            trueRes = trueRes.rxnsList;
                   
            testCase.verifyEqual(rxnsList, trueRes)
        end
        
        function testFixVariableNamesKinetics(testCase)
           
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            [xKinetic,strKinetic] = xlsread(xlsxFile, 'kinetics1');
            
            strKinetic = fixVariableNames(strKinetic, 'r', 'kinetics');
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResFixVariableNamesKinetics'));
            trueRes = trueRes.strKinetic;
                   
            testCase.verifyEqual(strKinetic, trueRes);
        end
        
        function testFixVariableNamesStoicCoeffs(testCase)
           
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_coef');
            [xKinetic,strKinetic] = xlsread(xlsxFile, 'kinetics1');
            
            strKinetic = fixVariableNames(strKinetic, 'r', 'kinetics');
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResFixVariableNamesStoicCoeffs'));
            trueRes = trueRes.strKinetic;
                   
            testCase.verifyEqual(strKinetic, trueRes);
        end
        
    end
end