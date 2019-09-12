classdef initialSamplerTest < matlab.unittest.TestCase
    
    methods (Test)
        function testInitialSampler1(testCase)
            seed = 1;
            rng(seed)
            
            % To generate the reaction files 
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.eigThreshold = 10^-5;
                        
            [isModelValid,model,strucIdx,xopt,tolScore,simFluxes] = initialSampler(ensemble);
            
            trueResModel = load(fullfile(currentPath{1}, 'testFiles', 'trueResModel_toy_model1.mat'));
            trueResModel = trueResModel.model;
            
            trueResIsModelValid = true;
            trueResModel.poolFactor = [];
            trueResStructIdx = 1;
            trueResXopt = 0;
            trueResTolScore = 0;
            trueResSimFluxes = 0;
           
            testCase.verifyEqual(trueResIsModelValid,isModelValid);
            testCase.verifyEqual(trueResModel,model);
            testCase.verifyEqual(trueResStructIdx,strucIdx);
            testCase.verifyEqual(trueResXopt,xopt);
            testCase.verifyEqual(trueResTolScore,tolScore);
            testCase.verifyEqual(trueResSimFluxes,simFluxes);
            
        end
        
        function testInitialSamplerWithRandom(testCase)
            seed = 1;
            rng(seed)
            
            % To generate the reaction files 
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');            
            xlsxFile = fullfile(currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.eigThreshold = 10^-5;
                        
            [isModelValid,model,strucIdx,xopt,tolScore,simFluxes] = initialSampler(ensemble);
           
            trueResModel = load(fullfile(currentPath{1}, 'testFiles', 'trueResModel_toy_model1_random.mat'));
            trueResModel = trueResModel.model;
            
            trueResIsModelValid = false;
            trueResModel.poolFactor = [];
            trueResStructIdx = 1;
            trueResXopt = 0;
            trueResTolScore = 0;
            trueResSimFluxes = 0;
           
            testCase.verifyEqual(trueResIsModelValid,isModelValid);
            testCase.verifyEqual(trueResModel,model);
            testCase.verifyEqual(trueResStructIdx,strucIdx);
            testCase.verifyEqual(trueResXopt,xopt);
            testCase.verifyEqual(trueResTolScore,tolScore);
            testCase.verifyEqual(trueResSimFluxes,simFluxes);
            
        end
    end
end