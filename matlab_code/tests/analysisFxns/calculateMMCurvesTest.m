classdef calculateMMCurvesTest < matlab.unittest.TestCase

    properties
        currentPath
        ensemble
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
        
        function loadEnsemble(testCase)
            testCase.ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'model100_OP2.mat'));
            testCase.ensemble = testCase.ensemble.ensemble;
        end
        
        function addReactionsToPath(testCase)
            addpath(fullfile(testCase.currentPath{1}, 'testFiles', 'reactions'));
        end
    end
    
   
    methods (Test)
        function testCalculateMMCurvesGTHOr(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [33];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_GTHOr_m_m_gthox_c.csv');
            resGthox = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_GTHOr_m_m_nadph_c.csv');
            resNadph = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_gthox.csv');
            trueResGthox = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_nadph.csv');
            trueResNadph = fileread(filepath);
            
            testCase.verifyEqual(trueResGthox, resGthox);
            testCase.verifyEqual(trueResNadph, resNadph);
        end
        
        function testCalculateMMCurvesGTHPi(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [32];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_GTHPi_m_m_gthrd_c.csv');
            resGthrd = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_GTHPi_m_m_h2o2_c.csv');
            resH2o2 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_gthrd.csv');
            trueResGthrd = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_h2o2.csv');
            trueResH2o2 = fileread(filepath);
            
            testCase.verifyEqual(trueResGthrd, resGthrd);
            testCase.verifyEqual(trueResH2o2, resH2o2);
        end
        
        function testCalculateMMCurvesG6PDH2NAD(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [11];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_G6PDH2_NAD_m_m_g6p_c.csv');
            resG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_G6PDH2_NAD_m_m_nad_c.csv');
            resNad = fileread(filepath);
            
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_g6p.csv');
            trueResG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_nad.csv');
            trueResNad = fileread(filepath);
            
            testCase.verifyEqual(trueResG6p, resG6p);
            testCase.verifyEqual(trueResNad, resNad);
        end
        
        function testCalculateMMCurvesG6PDH2NADP(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [12];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_G6PDH2_NADP_m_m_g6p_c.csv');
            resG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'r_R_G6PDH2_NADP_m_m_nadp_c.csv');
            resNadp = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_g6p.csv');
            trueResG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_nadp.csv');
            trueResNadp = fileread(filepath);
            
            testCase.verifyEqual(trueResG6p, resG6p);
            testCase.verifyEqual(trueResNadp, resNadp);
            
        end
       
    end
end
