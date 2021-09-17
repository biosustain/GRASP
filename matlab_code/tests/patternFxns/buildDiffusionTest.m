classdef buildDiffusionTest < matlab.unittest.TestCase
    
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
        function testBuildDiffusion1(testCase)
            
            reactionName = 'testDiffusion';
            strucIdx = 1;
            
            buildDiffusion(reactionName,strucIdx)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '1.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildDiffusion1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);            
        end
    end
end

