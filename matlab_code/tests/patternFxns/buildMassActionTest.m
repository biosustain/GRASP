classdef buildMassActionTest < matlab.unittest.TestCase
    
    properties
        currentPath
        tempReactionsFolder
    end
    
    methods(TestClassSetup)
        function createReactionsTempFolder(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
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
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildMassAction1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
        end
    end
end
