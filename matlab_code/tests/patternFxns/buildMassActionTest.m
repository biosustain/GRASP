classdef buildMassActionTest < matlab.unittest.TestCase
    
    methods (Test)
        function testBuildMassAction1(testCase)
            
            reactionName = 'testMassAction';
            strucIdx = 1;
            buildMassAction(reactionName,strucIdx)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildMassAction1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
        end
    end
end
