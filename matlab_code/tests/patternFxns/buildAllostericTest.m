classdef buildAllostericTest < matlab.unittest.TestCase

    methods (Test)
        function testBuildAllostericTestOneInhib(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9'};
            posEffectors = '';
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildAllostericOneInhib.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res);            
            
        end
        
        function testBuildAllostericTestOneInhibOneAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9'};
            posEffectors = {'m_m10'};
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildAllostericOneInhibOneAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res); 
        end
        
        function testBuildAllostericTestMultInhibOneAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9', 'm_m11'};
            posEffectors = {'m_m10'};
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildAllostericMultInhibOneAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res); 
        end
        
        function testBuildAllostericTestMultInhibMultAct(testCase)
            metList = {{'*A'}, {'*B'}, {'*C'}, {'*D'}, {'*P1'}, {'*P2'}, ...
                       {'*Q'}, {'*R'}};
            reactionName = 'r_r31';
            negEffectors = {'m_m9', 'm_m11'};
            posEffectors = {'m_m10'};
            buildAllosteric(metList,reactionName,negEffectors,posEffectors)
            
            
            currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
            tempReactionsFolder = fullfile(currentPath{1}, '..', '..', '..', 'temp', 'reactions');
            
            filepath = fullfile(tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(currentPath{1}, 'testFiles', 'trueResBuildAllostericMultInhibMultAct.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(trueRes,res); 
        end
    end
end

