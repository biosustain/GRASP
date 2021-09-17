classdef loadEnsembleStructureTest < matlab.unittest.TestCase

    properties
        currentPath
        relTol = 1e-2;
        absTol = 1e-4;        
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_allosteric_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random_1'),
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_freeExchange_1')};
            
            for i = 1:size(reactionsFolderList, 1)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
    end
       
    
    methods (Test)
        function testLoadEnsembleStructure1(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1.mat'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testLoadEnsembleStructure1WithAllostery(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_allosteric');
            ensemble = loadEnsembleStructure(xlsxFile);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_allosteric.mat'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testLoadEnsembleStructure1WithRandom(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random');
            ensemble = loadEnsembleStructure(xlsxFile);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_random.mat'));            
            trueRes = trueRes.ensemble;
                                              
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
        
        function testLoadEnsembleStructureWithFreeExchange(testCase)
            % This test is expected to fail
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_freeExchange');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model1_freeExchange.mat'));
            trueRes = trueRes.ensemble;            
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));         
        end
       
        function testLoadEnsembleStructureInputValidationGeneral(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_0');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationRxns(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_1');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationRxnsRows(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_1_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationMets(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_2');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationMetsRows(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_2_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationMeasRates(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_3');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationThermoRxns(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_4');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end

        function testLoadEnsembleStructureInputValidationThermoRxnsRows(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_4_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationThermoMets(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_5');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationThermoMetsRows(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_5_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationMetsData(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_6');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end

        function testLoadEnsembleStructureInputValidationMetsDataRow(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_6_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationProtData(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_7');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationProtDataRows(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_7_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationKinetics(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_8');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationKineticsRow(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_input_valid_8_row');
            testCase.verifyError(@()loadEnsembleStructure(xlsxFile), '');            
            
        end
        
        function testLoadEnsembleStructureInputValidationMetsActive(testCase)
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model4');
            ensemble = loadEnsembleStructure(xlsxFile);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResLoadEnsemble_toy_model4.mat'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end

