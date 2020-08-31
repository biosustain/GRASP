classdef buildOdeFxnTest < matlab.unittest.TestCase
            
    properties
        currentPath
        reactionsFolder
    end
    
    methods(TestClassSetup)
        function createReactionsFolder(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
            testCase.reactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1');
            mkdir(testCase.reactionsFolder);
            
            kineticsFile = fullfile(testCase.currentPath{1}, 'testFiles','toy_model1_random2_Kinetics1.m');
            copyfile(kineticsFile, testCase.reactionsFolder);
        end
    end
 
    methods(TestClassTeardown)
        function removeReactionsFolder(testCase)           
            if exist(testCase.reactionsFolder, 'dir')
                rmdir(testCase.reactionsFolder, 's');
            end
        end
    end
    
    
    methods (Test)
        function testBuildOdeFxn1(testCase)
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1_random2.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            kineticFxn = 'toy_model1_random2_Kinetics1';
            strucIdx = 1;
            
            buildOdeFxn(ensemble,kineticFxn,strucIdx);

            res = fileread(fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', [ensemble.description, '_', num2str(strucIdx)], [kineticFxn,'_ode.m']));
            trueRes = fileread(fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_toy_model1_random2_Kinetics1_ode.m'));
            
            testCase.verifyEqual(res, trueRes);                   
        end
    end
end

