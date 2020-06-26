classdef buildEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_allosteric2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_new_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_no_promiscuous2_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_debug_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model5_dGs_1')};
            
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
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');
            trueRes.LPSolver = 'gurobi';
                   
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildEnsembleRandom'));
            trueRes = trueRes.ensemble;
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');
            trueRes.LPSolver = 'linprog';
                   
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');   
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes, 'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');   
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');       
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes, 'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');       
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');       
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            
            trueRes = rmfield(trueRes,'rxnMetLinks');
            trueRes = rmfield(trueRes,'freeVars');
            ensemble = rmfield(ensemble,'freeVars');       
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
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
            trueRes.LPSolver = 'gurobi';
            
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(ensemble, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
	end
end




