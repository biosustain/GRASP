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
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHOr_m_m_gthox_c.csv');
            resGthox = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHOr_m_m_nadph_c.csv');
            resNadph = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_gthox.csv');
            trueResGthox = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_nadph.csv');
            trueResNadph = fileread(filepath);
            
            testCase.verifyEqual(resGthox, trueResGthox);
            testCase.verifyEqual(resNadph, trueResNadph);
        end
        
        function testCalculateMMCurvesGTHPi(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [32];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHPi_m_m_gthrd_c.csv');
            resGthrd = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHPi_m_m_h2o2_c.csv');
            resH2o2 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_gthrd.csv');
            trueResGthrd = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_h2o2.csv');
            trueResH2o2 = fileread(filepath);
            
            testCase.verifyEqual(resGthrd, trueResGthrd);
            testCase.verifyEqual(resH2o2, trueResH2o2);
        end
        
        function testCalculateMMCurvesG6PDH2NAD(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [11];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NAD_m_m_g6p_c.csv');
            resG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NAD_m_m_nad_c.csv');
            resNad = fileread(filepath);
            
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_g6p.csv');
            trueResG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_nad.csv');
            trueResNad = fileread(filepath);
            
            testCase.verifyEqual(resG6p, trueResG6p);
            testCase.verifyEqual(resNad, trueResNad);
        end
        
        function testCalculateMMCurvesG6PDH2NADP(testCase)
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 100;
            structIdx = 1;
            saturatingConc = 10^2*10^6;
            substrateRange = logspace(-15,6);
            rxnList = [12];
            calculateMMCurves(outputFolder, testCase.ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_g6p_c.csv');
            resG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_nadp_c.csv');
            resNadp = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_6pgl_c.csv');
            res6PGL = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_nadph_c.csv');
            resNadph = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_g6p.csv');
            trueResG6p = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_nadp.csv');
            trueResNadp = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_6pgl.csv');
            trueRes6PGL = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_nadph.csv');
            trueResNadph = fileread(filepath);
            
            testCase.verifyEqual(resG6p, trueResG6p);
            testCase.verifyEqual(resNadp, trueResNadp);
            testCase.verifyEqual(res6PGL, trueRes6PGL);
            testCase.verifyEqual(resNadph, trueResNadph);
            
        end
        
        function testCalculateMMCurvesInhibitors(testCase)
            
            ensembleLocal = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model.mat'));
            ensembleLocal = ensembleLocal.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 5;
            structIdx = 1;
            saturatingConc = 10^4;
            substrateRange = logspace(-6, 4);
                        
            calculateMMCurves(outputFolder, ensembleLocal, numModels, structIdx, saturatingConc, substrateRange);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r1_m_m3.csv');
            resR1M3 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r1_m_m6.csv');
            resR1M6 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r3_m_m1.csv');
            resR3M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r3_m_m7.csv');
            resR3M7 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m3.csv');
            trueResR1M3 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m6.csv');
            trueResR1M6 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m1.csv');
            trueResR3M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m7.csv');
            trueResR3M7 = fileread(filepath);
            
            testCase.verifyEqual(resR1M3, trueResR1M3);
            testCase.verifyEqual(resR1M6, trueResR1M6);
            testCase.verifyEqual(resR3M1, trueResR3M1);
            testCase.verifyEqual(resR3M7, trueResR3M7);
        end
        
        function testCalculateMMCurvesWater(testCase)
            
            ensembleLocal = load(fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_water.mat'));
            ensembleLocal = ensembleLocal.ensemble;
            
            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles');
            numModels = 5;
            structIdx = 1;
            
            saturatingConc = 10^4;
            substrateRange = logspace(-6, 4);
            rxnList = [1, 3, 6];

            calculateMMCurves(outputFolder, ensembleLocal, numModels, structIdx, saturatingConc, substrateRange, rxnList);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r1_m_m3.csv');
            resR1M3 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r1_m_m6.csv');
            resR1M6 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r3_m_m1.csv');
            resR3M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r3_m_m7.csv');
            resR3M7 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r6_m_m1.csv');
            resR6M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r6_m_m10.csv');
            resR6M10 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m3_water.csv');
            trueResR1M3 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m6_water.csv');
            trueResR1M6 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m1_water.csv');
            trueResR3M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m7_water.csv');
            trueResR3M7 = fileread(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r6_m_m1_water.csv');
            trueResR6M1 = fileread(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r6_m_m10_water.csv');
            trueResR6M10 = fileread(filepath);
            
            testCase.verifyEqual(resR1M3, trueResR1M3);
            testCase.verifyEqual(resR1M6, trueResR1M6);
            testCase.verifyEqual(resR3M1, trueResR3M1);
            testCase.verifyEqual(resR3M7, trueResR3M7);
            testCase.verifyEqual(resR6M1, trueResR6M1);
            testCase.verifyEqual(resR6M10, trueResR6M10);
            testCase.verifyTrue(~exist(fullfile(testCase.currentPath{1}, 'testFiles', 'r_r1_m_m0.csv'), 'file'));
            testCase.verifyTrue(~exist(fullfile(testCase.currentPath{1}, 'testFiles', 'r_r1_m_m0.csv'), 'file'));
        end
       
    end
end
