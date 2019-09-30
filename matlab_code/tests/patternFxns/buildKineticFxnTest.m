classdef buildKineticFxnTest < matlab.unittest.TestCase
        
    properties
        currentPath
        reactionsFolder
    end
    
    methods(TestClassSetup)
        function createReactionsFolder(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
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
            
            [rxnMetLinks,freeVars,metsActive] = buildKineticFxn(ensemble,kineticFxn,strucIdx);
            
            trueResRxnMetLinks = {'r_r1', 'r_r2', 'r_r3', 'r_r4', 'r_r5', 'r_r6', ...
                                  {'m_m6', 'r_r7'}, {'m_m6', 'r_r8'}, {'m_m7', 'r_r9'}, ...
                                  {'m_m5', 'r_r10'}, 'm_m8r_r11', {'m_m9', 'r_r12'}, ...
                                  []};
                              
            trueResFreeVars = {'m_m5'; 'm_m6'; 'm_m7'; 'm_m8'; 'm_m9';
                               'm_m10'; 'm_m11'; 'r_r1'; 'r_r2'; 'r_r3';
                               'r_r4'; 'r_r5'; 'r_r6'; 'r_r7'; 'r_r8';
                               'r_r9'; 'r_r10'; 'r_r11'; 'r_r12'};
                           
            trueResMetsActive = [5; 6; 7; 8; 9; 10; 11];
            
            
            testCase.verifyEqual(trueResRxnMetLinks,rxnMetLinks);            
            testCase.verifyEqual(trueResFreeVars,freeVars);            
            testCase.verifyEqual(trueResMetsActive, metsActive);            
        end
    end
end

