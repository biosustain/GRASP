classdef buildFreeExchangeTest < matlab.unittest.TestCase
    
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
        function testBuildFreeExchange1(testCase)
            
            reactionName = 'testFreeExchange';
            strucIdx = 1;
            buildFreeExchange(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildFreeExchange1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);            
        end
    end
end
