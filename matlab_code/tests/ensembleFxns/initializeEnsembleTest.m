classdef initializeEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
    
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolderList = {fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model3_isoenzymes_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model4_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model4_inhibitors_1'), ...
                                   fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_nick_1')};
            
            for i = 1:size(reactionsFolderList, 2)
                disp( exist(reactionsFolderList{i}, 'dir'));
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
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleWithRandom(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1_random.mat');
            ensemble = load(filepath);
            ensemble = ensemble.ensemble;
            popIdx = 1;
            verbose = 1;
            
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_random.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleIsoenzymes(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model3_isoenzymes');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model3_isoenzymes.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
         
        function testInitializeEnsembleToyModel4(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model4');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model4.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleToyModel4Inhibitors(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model4_inhibitors');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model4_inhibitors.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
        
        function testInitializeEnsembleToyModel1Nick(testCase)
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_nick');
            popIdx = 1;
            verbose = 1;
            
            ensemble = loadEnsembleStructure(filepath);   
            ensemble = initializeEnsemble(ensemble,popIdx,verbose);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResInitializedEnsemble_toy_model1_nick.mat'));
            testCase.verifyEqual(trueRes.ensemble,ensemble);   
            
        end
    end
end

