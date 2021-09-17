classdef initialSamplerTest < matlab.unittest.TestCase

    properties
        currentPath
        relTol = 1e-4;
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
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random_1')};
            
            for i = 1:size(reactionsFolderList, 1)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
    end
       
    
    methods (Test)
        function testInitialSampler1(testCase)
            seed = 1;
            rng(seed)
            
            % To generate the reaction files 
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);
         
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.sampler = 'GRASP';
            ensemble.eigThreshold = 10^-5;
            modelI = 1;
            [~,ensemble.metsLI] = rref(ensemble.Sred');
                        
            [isModelValid,model,strucIdx,xopt,tolScore,simFluxes] = initialSampler(ensemble, modelI);

            trueResModel = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResModel_toy_model1.mat'));
            trueResModel = trueResModel.model;
            
            trueResIsModelValid = true;
            trueResModel.poolFactor = [];
            trueResStructIdx = 1;
            trueResXopt = 0;
            trueResTolScore = 0;
            trueResSimFluxes = 0;
                       
            testCase.verifyEqual(isModelValid, trueResIsModelValid);
            testCase.verifyEqual(strucIdx, trueResStructIdx);
            testCase.verifyEqual(xopt, trueResXopt);
            testCase.verifyEqual(tolScore, trueResTolScore);
            testCase.verifyEqual(simFluxes, trueResSimFluxes);
            testCase.verifyThat(model, matlab.unittest.constraints.IsEqualTo(trueResModel, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
        
        function testInitialSamplerWithRandom(testCase)
            seed = 1;
            rng(seed)
            
            % To generate the reaction files 
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            ensemble = loadEnsembleStructure(xlsxFile);

            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.sampler = 'GRASP';
            ensemble.eigThreshold = 10^-5;
            ensemble.freeVars{end+1} = 'r_r13';
            ensemble.LPSolver = 'gurobi';
            modelI = 1;
            [~,ensemble.metsLI] = rref(ensemble.Sred');
                        
            [isModelValid,model,strucIdx,xopt,tolScore,simFluxes] = initialSampler(ensemble, modelI);

            trueResModel = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResModel_toy_model1_random.mat'));
            trueResModel = trueResModel.model;
            
            trueResIsModelValid = true;
            trueResModel.poolFactor = [];
            trueResStructIdx = 1;
            trueResXopt = 0;
            trueResTolScore = 0;
            trueResSimFluxes = 0;           

            testCase.verifyEqual(isModelValid, trueResIsModelValid);
            testCase.verifyEqual(strucIdx, trueResStructIdx);
            testCase.verifyEqual(xopt, trueResXopt);
            testCase.verifyEqual(tolScore, trueResTolScore);
            testCase.verifyEqual(simFluxes, trueResSimFluxes);
            testCase.verifyThat(model, matlab.unittest.constraints.IsEqualTo(trueResModel, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
        
        function testInitialSamplerWithConstEffector(testCase)
            seed = 1;
            rng(seed)
            
            % To generate the reaction files 
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_const_effector');
            ensemble = loadEnsembleStructure(xlsxFile);
            ensemble = initializeEnsemble(ensemble,1,1);
            ensemble.eigThreshold = 10^-5;
                        
            addKineticFxnsToPath(ensemble);
            maxNumberOfSamples = 100;

            ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples);
            
            modelI = 1;
            [isModelValid,model,strucIdx,xopt,tolScore,simFluxes] = initialSampler(ensemble, modelI);

            trueResModel = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResModel_toy_model1_const_effector.mat'));
            trueResModel = trueResModel.model;
                        
            trueResIsModelValid = true;
            trueResModel.poolFactor = [];
            trueResStructIdx = 1;
            trueResXopt = 0;
            trueResTolScore = 0;
            trueResSimFluxes = 0;
           
            testCase.verifyEqual(isModelValid, trueResIsModelValid);
            testCase.verifyEqual(strucIdx, trueResStructIdx);
            testCase.verifyEqual(xopt, trueResXopt);
            testCase.verifyEqual(tolScore, trueResTolScore);
            testCase.verifyEqual(simFluxes, trueResSimFluxes);
            testCase.verifyThat(model, matlab.unittest.constraints.IsEqualTo(trueResModel, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
    end
end