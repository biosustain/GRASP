classdef buildEnsembleTest < matlab.unittest.TestCase

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

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model_numbers_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_linprog_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_allosteric2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_new_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_no_promiscuous2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_debug_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_dGs_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_parallel_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'methionine_cycle_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'methionine_cycle_fmincon_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'methionine_cycle_atp_change_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'methionine_cycle_pedro_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_debug_1')};
            
            for i = 1:size(reactionsFolderList, 2)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
        
        function removeOutputFolder(testCase)           

            outputFolder = fullfile(testCase.currentPath{1}, 'testFiles', 'testDir1');
            
            if exist(outputFolder, 'dir')
                rmdir(outputFolder, 's');
            end
        end
    end
    
    
    methods (Test)
        function testBuildEnsembleNumbers(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_numbers';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleNumbers'));
            trueRes = trueRes.ensemble;
                               
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleRandom(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_random2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleRandom'));
            trueRes = trueRes.ensemble;
                               
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end

        function testBuildEnsembleRandomLinprog(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_random2_linprog';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleRandomLinprog'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within',matlab.unittest.constraints.RelativeTolerance(1.7e-1) | matlab.unittest.constraints.AbsoluteTolerance(1.4e-1)));
        end
               
        function testBuildEnsembleAllosteric(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_allosteric2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleAllosteric'));
            trueRes = trueRes.ensemble;
            trueRes.sampler = 'GRASP';  
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleAllStable(testCase)
            % All models are expected to be stable
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_new';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
 
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleAllStable'));
            trueRes = trueRes.ensemble;          
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleNoPromiscuous(testCase)
            % 10-20% of all models are expected to be unstable
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_no_promiscuous2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleNoPromiscuous'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleLargeModel(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model5_debug';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleLargeModel'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(1e-1)));
        end
        
        function testBuildEnsembleLargeModeldGs(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model5_dGs';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleLargeModeldGs'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
        
        function testBuildEnsembleLargeModelParallel(testCase)
            
           
            modelID = 'toy_model5_parallel';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            testing = true;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold, testing);
          
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleLargeModelParallel'));
            trueRes = trueRes.ensemble;

            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(1e-1)));
        end
        
        function testBuildEnsembleCreateDir(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'toy_model1_random2';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', 'testDir1', 'testDir2', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensembleTemp = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            ensemble = load(outputFile);
            ensemble = ensemble.ensemble;
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleCreateDir'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end

        function testBuildEnsembleExample(testCase)
            
            seed = 1;
            rng(seed)

            
            modelID = 'toy_model';
            inputFile = fullfile(testCase.currentPath{1}, '..', '..', '..', 'io', 'input', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', 'testDir1', 'testDir2', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensembleTemp = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
           
            ensemble = load(outputFile);
            ensemble = ensemble.ensemble;

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleExample'));
            trueRes = trueRes.ensemble;  
                        
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleExampleABCNlopt(testCase)
            
            seed = 1;
            rng(seed)

            
            modelID = 'methionine_cycle';
            inputFile = fullfile(testCase.currentPath{1}, '..', '..', '..', 'io', 'input', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', 'testDir1', 'testDir2', [modelID, '.mat']);
            
            maxNumberOfSamples = 100;
            eigThreshold = 10^-5;
            
            ensemble= buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleExampleABCNLopt'));
            trueRes = trueRes.ensemble;
                        
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
        end
        
        function testBuildEnsembleExampleABCFmincon(testCase)
            
            seed = 1;
            rng(seed)

            
            modelID = 'methionine_cycle_fmincon';
            inputFile = fullfile(testCase.currentPath{1}, '..', '..', '..', 'io', 'input', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', 'testDir1', 'testDir2', [modelID, '.mat']);
            
            maxNumberOfSamples = 100;
            eigThreshold = 10^-5;
            
            ensembleTemp = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);
            
            ensemble = load(outputFile);
            ensemble = ensemble.ensemble;

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleExampleABCFmincon'));
            trueRes = trueRes.ensemble;
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1.7e-2) | matlab.unittest.constraints.AbsoluteTolerance(2e-3)));
        end
        
        function testBuildEnsembleMethionineATPchange(testCase)
            
            seed = 1;
            rng(seed)
            
            modelID = 'methionine_cycle_atp_change';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResTestBuildEnsembleMethionineATPchange'));
            trueRes = trueRes.ensemble;
       
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
        
        function testBuildEnsembleMethioninePedroParallel(testCase)
                        
            modelID = 'methionine_cycle_pedro';
            inputFile = fullfile(testCase.currentPath{1}, 'testFiles', modelID);
            outputFile = fullfile(testCase.currentPath{1}, 'testFiles', [modelID, '.mat']);
            
            maxNumberOfSamples = 1000;
            eigThreshold = 10^-5;
            testing = true;
            
            ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold, testing);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResTestBuildEnsembleMethioninePedroParallel'));
            trueRes = trueRes.ensemble;
       
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(1.5e-3)));
        end
	end
end




