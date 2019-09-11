classdef buildFreeExchangeTest < matlab.unittest.TestCase
    
    methods (Test)
        function testBuildFreeExchange1(testCase)
            
            reactionName = 'testFreeExchange';
            strucIdx = 1;
            buildFreeExchange(reactionName,strucIdx)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildFreeExchange1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
        end
    end
end
