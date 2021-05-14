classdef calculateMMCurvesTest < matlab.unittest.TestCase

    properties
        currentPath
        ensemble
        relTol = 1e-10;
        absTol = 1e-10;
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
            resGthox = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHOr_m_m_nadph_c.csv');
            resNadph = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_gthox.csv');
            trueResGthox = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHOr_nadph.csv');
            trueResNadph = readmatrix(filepath);
            
            testCase.verifyThat(resGthox, matlab.unittest.constraints.IsEqualTo(trueResGthox, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resNadph, matlab.unittest.constraints.IsEqualTo(trueResNadph, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
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
            resGthrd = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_GTHPi_m_m_h2o2_c.csv');
            resH2o2 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_gthrd.csv');
            trueResGthrd = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_GTHPi_h2o2.csv');
            trueResH2o2 = readmatrix(filepath);
            
            testCase.verifyThat(resGthrd, matlab.unittest.constraints.IsEqualTo(trueResGthrd, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resH2o2, matlab.unittest.constraints.IsEqualTo(trueResH2o2, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
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
            resG6p = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NAD_m_m_nad_c.csv');
            resNad = readmatrix(filepath);
            
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_g6p.csv');
            trueResG6p = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NAD_nad.csv');
            trueResNad = readmatrix(filepath);
            
            testCase.verifyThat(resG6p, matlab.unittest.constraints.IsEqualTo(trueResG6p, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resNad, matlab.unittest.constraints.IsEqualTo(trueResNad, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
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
            resG6p = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_nadp_c.csv');
            resNadp = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_6pgl_c.csv');
            res6PGL = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'putida_OP2_r_R_G6PDH2_NADP_m_m_nadph_c.csv');
            resNadph = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_g6p.csv');
            trueResG6p = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_nadp.csv');
            trueResNadp = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_6pgl.csv');
            trueRes6PGL = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'OP2_G6PDH2_NADP_nadph.csv');
            trueResNadph = readmatrix(filepath);

            testCase.verifyThat(resG6p, matlab.unittest.constraints.IsEqualTo(trueResG6p, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resNadp, matlab.unittest.constraints.IsEqualTo(trueResNadp, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(res6PGL, matlab.unittest.constraints.IsEqualTo(trueRes6PGL, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resNadph, matlab.unittest.constraints.IsEqualTo(trueResNadph, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
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
            resR1M3 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r1_m_m6.csv');
            resR1M6 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r3_m_m1.csv');
            resR3M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model_r_r3_m_m7.csv');
            resR3M7 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m3.csv');
            trueResR1M3 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m6.csv');
            trueResR1M6 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m1.csv');
            trueResR3M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m7.csv');
            trueResR3M7 = readmatrix(filepath);

            testCase.verifyThat(resR1M3, matlab.unittest.constraints.IsEqualTo(trueResR1M3, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR1M6, matlab.unittest.constraints.IsEqualTo(trueResR1M6, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR3M1, matlab.unittest.constraints.IsEqualTo(trueResR3M1, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR3M7, matlab.unittest.constraints.IsEqualTo(trueResR3M7, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
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
            resR1M3 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r1_m_m6.csv');
            resR1M6 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r3_m_m1.csv');
            resR3M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r3_m_m7.csv');
            resR3M7 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r6_m_m1.csv');
            resR6M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2_r_r6_m_m10.csv');
            resR6M10 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m3_water.csv');
            trueResR1M3 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r1_m_m6_water.csv');
            trueResR1M6 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m1_water.csv');
            trueResR3M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r3_m_m7_water.csv');
            trueResR3M7 = readmatrix(filepath);
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r6_m_m1_water.csv');
            trueResR6M1 = readmatrix(filepath);
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_r_r6_m_m10_water.csv');
            trueResR6M10 = readmatrix(filepath);
            
            testCase.verifyThat(resR1M3, matlab.unittest.constraints.IsEqualTo(trueResR1M3, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR1M6, matlab.unittest.constraints.IsEqualTo(trueResR1M6, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR3M1, matlab.unittest.constraints.IsEqualTo(trueResR3M1, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR3M7, matlab.unittest.constraints.IsEqualTo(trueResR3M7, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR6M1, matlab.unittest.constraints.IsEqualTo(trueResR6M1, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            testCase.verifyThat(resR6M10, matlab.unittest.constraints.IsEqualTo(trueResR6M10, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
            testCase.verifyTrue(~exist(fullfile(testCase.currentPath{1}, 'testFiles', 'r_r1_m_m0.csv'), 'file'));
            testCase.verifyTrue(~exist(fullfile(testCase.currentPath{1}, 'testFiles', 'r_r1_m_m0.csv'), 'file'));
        end
       
    end
end
