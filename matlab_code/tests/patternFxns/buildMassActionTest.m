classdef buildMassActionTest < matlab.unittest.TestCase
    
    properties
        currentPath
        tempReactionsFolder
    end
    
    methods(TestClassSetup)
        function createReactionsTempFolder(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
            testCase.tempReactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            mkdir(testCase.tempReactionsFolder);
        end
    end
 
    methods(TestClassTeardown)
        function removeReactionsTempFolder(testCase)           
            if exist(testCase.tempReactionsFolder, 'dir')
                rmdir(testCase.tempReactionsFolder, 's');
            end
        end
    end
 
    
    methods (Test)
        function testBuildMassAction1(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = fileread(filepath);
            
            addpath(testCase.tempReactionsFolder);
            SC = [1];
            S = [0.1];
            PC = [1];
            P = [0.2];
            K = [1, 2];
            v = testMassAction1(SC,S,PC,P,K);
            
            trueResV = -0.3;           
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.m');
            trueRes = fileread(filepath);
            
            testCase.verifyEqual(res, trueRes); 
            testCase.verifyThat(v, matlab.unittest.constraints.IsEqualTo(trueResV, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
        
        function testBuildMassActionOrder(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = fileread(filepath);
            
            addpath(testCase.tempReactionsFolder);
            SC = [1; 2];
            S = [0.1; 2];
            PC = [1];
            P = [0.2];
            K = [1, 2];
            v = testMassAction1(SC,S,PC,P,K);
            
            trueResV = 0;           
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.m');
            trueRes = fileread(filepath);
            
            testCase.verifyEqual(res, trueRes); 
            testCase.verifyThat(v , matlab.unittest.constraints.IsEqualTo(trueResV, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
            
        end

        function testBuildMassActionOrderCoef(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = fileread(filepath);
            
            addpath(testCase.tempReactionsFolder);
            SC = [1; 2];
            S = [3; 2];
            PC = [2; 1];
            P = [0.2; 1.5];
            K = [1, 2];
            v = testMassAction1(SC,S,PC,P,K);
            
            trueResV = 11.8800;           
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.m');
            trueRes = fileread(filepath);
            
            testCase.verifyEqual(res, trueRes); 
            testCase.verifyThat(v, matlab.unittest.constraints.IsEqualTo(trueResV, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
        
        function testBuildMassActionControlAnalysis(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = fileread(filepath);
            
            addpath(testCase.tempReactionsFolder);
            load(fullfile(testCase.currentPath{1}, 'testFiles', 'testBuilMassActionControlAnalysis.mat'));
            v = testMassAction1(SC,S,PC,P,K);

            trueResV = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResVBuilMassActionControlAnalysis.mat'));
            trueResV = trueResV.trueResV;
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.m');
            trueRes = fileread(filepath);
            
            testCase.verifyEqual(res, trueRes); 
            testCase.verifyThat(v, matlab.unittest.constraints.IsEqualTo(trueResV, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
        
        function testBuildMassActionControlAnalysis2(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = fileread(filepath);
            
            addpath(testCase.tempReactionsFolder);
            load(fullfile(testCase.currentPath{1}, 'testFiles', 'testBuilMassActionControlAnalysis.mat'));
            P = [P; 2*ones(1,7)];
            PC = [1*ones(1,7); 2*ones(1,7)];
            v = testMassAction1(SC,S,PC,P,K);
            
            trueResV = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResVBuilMassActionControlAnalysis2.mat'));
            trueResV = trueResV.trueResV;
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.m');
            trueRes = fileread(filepath);
            
            testCase.verifyEqual(res, trueRes); 
            testCase.verifyThat(v, matlab.unittest.constraints.IsEqualTo(trueResV, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
    end
end
