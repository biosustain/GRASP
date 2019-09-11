classdef buildDiffusionTest < matlab.unittest.TestCase
    
    methods (Test)
        function testBuildDiffusion1(testCase)
            
            reactionName = 'testDiffusion';
            strucIdx = 1;
            buildDiffusion(reactionName,strucIdx)
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildDiffusion1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
        end
    end
end

