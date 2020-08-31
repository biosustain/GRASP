classdef buildAllostericTest < matlab.unittest.TestCase
    
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
        function testBuildAllostericTestOneInhib(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9'};
            posEffectors = '';
            
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildAllostericOneInhib.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);            
            
        end
        
        function testBuildAllostericTestOneInhibOneAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9'};
            posEffectors = {'m_m10'};
            
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildAllostericOneInhibOneAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res,trueRes); 
        end
        
        function testBuildAllostericTestMultInhibOneAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9', 'm_m11'};
            posEffectors = {'m_m10'};
            
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildAllostericMultInhibOneAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes); 
        end
        
        function testBuildAllostericTestMultInhibMultAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9', 'm_m11'};
            posEffectors = {'m_m10'};
            
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildAllostericMultInhibMultAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes); 
        end
    end
end

