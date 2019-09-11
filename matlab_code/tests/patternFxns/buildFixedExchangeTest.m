classdef buildFixedExchangeTest < matlab.unittest.TestCase
    
    methods (Test)
        function testBuildFixedExchange1(testCase)
            
            reactionName = 'testFixedExchange';
            strucIdx = 1;
            buildFixedExchange(reactionName,strucIdx)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildFixedExchange1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
        end
    end
end
