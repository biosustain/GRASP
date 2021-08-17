classdef initializeEnsembleTest < matlab.unittest.TestCase

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

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model3_isoenzymes_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model4_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model4_inhibitors_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_nick_1')};
            
            for i = 1:size(reactionsFolderList, 2)
                if exist(reactionsFolderList{i}, 'dir')
                    rmdir(reactionsFolderList{i}, 's');
                end
            end
        end
    end    
 
    methods (Test)
        function testInitializeEnsemble(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.sampler = 'GRASP';
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1.mat'));
            
            trueRes = trueRes.ensemble;
            trueRes.sampler = 'GRASP';
          
             testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
        
        function testInitializeEnsembleWithRandom(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            ensemble.sampler = 'GRASP';
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_random.mat'));
            
            trueRes = trueRes.ensemble;
            trueRes.sampler = 'GRASP'; 
            
             testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
            
        end
        
        function testInitializeEnsembleIsoenzymes(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model3_isoenzymes');
            popIdx = 1;
            verbose = 1;
                        
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model3_isoenzymes.mat'));
            trueRes = trueRes.ensemble;   
       
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
            
        end
         
        function testInitializeEnsembleToyModel4(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model4');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model4.mat'));
            trueRes = trueRes.ensemble; 
            
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
            
        end
        
        function testInitializeEnsembleToyModel4Inhibitors(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model4_inhibitors');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model4_inhibitors.mat'));
            trueRes = trueRes.ensemble;
                       
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
            
        end
        
        function testInitializeEnsembleToyModel1Nick(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_nick');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_nick.mat'));
            trueRes = trueRes.ensemble;       
                       
            testCase.verifyThat(ensemble, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(8e-3) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));
            
        end
    end
end

