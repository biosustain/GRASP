classdef buildKineticFxnTest < matlab.unittest.TestCase
        
    properties
        currentPath
        reactionsFolder
    end
    
    methods(TestClassSetup)
        function createReactionsFolder(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
            testCase.reactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_1');
            mkdir(testCase.reactionsFolder);
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
        function testKineticFxn1(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            kineticFxn = 'toy_model1_Kinetics1';
            strucIdx = 1;
            
            freeVars = buildKineticFxn(ensemble,kineticFxn,strucIdx);
            
            kineticFxn = fileread(fullfile(testCase.reactionsFolder, 'toy_model1_Kinetics1.m'));
            
            trueResKineticFxn = fileread(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildKineticFxn1.m'));
                                        
            trueResFreeVars = {'m_m5'; 'm_m6'; 'm_m7'; 'm_m8'; 'm_m9';
                               'm_m10'; 'm_m11'; 'r_r1'; 'r_r2'; 'r_r3';
                               'r_r4'; 'r_r5'; 'r_r6'; 'r_r7'; 'r_r8';
                               'r_r9'; 'r_r10'; 'r_r11'; 'r_r12'; 'r_r13'};
                                       
            testCase.verifyEqual(kineticFxn, trueResKineticFxn);  
            testCase.verifyEqual(freeVars, trueResFreeVars);            
        end
        
        function testKineticFxnStoicCoef(testCase)
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model_stoic_coeff.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            kineticFxn = 'toy_model1_Kinetics1';
            strucIdx = 1;
            
            freeVars = buildKineticFxn(ensemble,kineticFxn,strucIdx);
            
            kineticFxn = fileread(fullfile(testCase.reactionsFolder, 'toy_model1_Kinetics1.m'));
            
            trueResKineticFxn = fileread(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildKineticFxnStoicCoef.m'));
            
        
            trueResFreeVars = {'m_m5'; 'm_m6'; 'm_m7'; 'm_m8'; 'm_m9';
                               'm_m10'; 'm_m11'; 'r_r1'; 'r_r2'; 'r_r3';
                               'r_r4'; 'r_r5'; 'r_r6'; 'r_r7'; 'r_r8';
                               'r_r9'; 'r_r10'; 'r_r11'; 'r_r12'; 'r_r13'};
                                                 
            testCase.verifyEqual(kineticFxn, trueResKineticFxn);  
            testCase.verifyEqual(freeVars, trueResFreeVars);            
        end
    end
end

