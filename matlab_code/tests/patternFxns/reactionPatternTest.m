classdef reactionPatternTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testReactionPattern1(testCase)
            
            patternName = ['substrateInhibOrderedBiBi'];
            reactionName = 'r_r1';
            flag = 1;
            strucIdx = 1;
            promiscuousRxnI = 0;
            
            
            [revMatrix,forwardFlux,metList] = reactionPattern(patternName,reactionName,flag,strucIdx, promiscuousRxnI);
            
            trueResRevMatrix = [1,0,1,1,1,1];
            trueResForwardFlux = [1,2;1,3;2,4;4,5;5,6;6,1];
            trueResMetList = {'*A','*B','*I','*P','*Q'};
           
            testCase.verifyEqual(trueResRevMatrix, revMatrix);  
            testCase.verifyEqual(trueResForwardFlux, forwardFlux);  
            testCase.verifyEqual(trueResMetList, metList);  
            
        end
        
        function testReactionPatternOrdereBiBi(testCase)
            
            patternName = ['orderedBiBi'];
            reactionName = 'r_r4';
            flag = 1;
            strucIdx = 1;
            promiscuousRxnI = 0;
            
            
            [revMatrix,forwardFlux,metList] = reactionPattern(patternName,reactionName,flag,strucIdx, promiscuousRxnI);
            
            trueResRevMatrix = [1,1,1,1,1];
            trueResForwardFlux = [1,2;2,3;3,4;4,5;5,1];
            trueResMetList = {'*A','*B','*P','*Q'};
           
            testCase.verifyEqual(trueResRevMatrix, revMatrix);  
            testCase.verifyEqual(trueResForwardFlux, forwardFlux);  
            testCase.verifyEqual(trueResMetList, metList);  
            
        end
        
        function testReactionPatternPromiscuous(testCase)
            
            patternName = ['OrdPromiscCompInhibIndep'];
            reactionName = 'r_r3';
            flag = 1;
            strucIdx = 1;
            promiscuousRxnI = 1;
            
            
            [revMatrix,forwardFlux,metList] = reactionPattern(patternName,reactionName,flag,strucIdx, promiscuousRxnI);
            
            trueResRevMatrix = [0,0,0,0,0,0,1,0,1,1,1,1;1,0,1,1,1,1,0,0,0,0,0,0];
            trueResForwardFlux = [1,2;2,3;2,4;4,5;5,6;6,1;1,7;7,8;7,9;9,10;10,11;11,1];
            trueResMetList = {'*A','*B','*C','*D','*I','*P1','*P2','*Q','*R'};
           
            testCase.verifyEqual(trueResRevMatrix, revMatrix);  
            testCase.verifyEqual(trueResForwardFlux, forwardFlux);  
            testCase.verifyEqual(trueResMetList, metList);  
            
        end
        
        function testReactionPatternRandom(testCase)
            
            patternName = ['randomBiBi'];
            reactionName = 'r_r4';
            flag = 1;
            strucIdx = 1;
            promiscuousRxnI = 0;
            
            
            [revMatrix,forwardFlux,metList] = reactionPattern(patternName,reactionName,flag,strucIdx, promiscuousRxnI);
            
            trueResRevMatrix = [0,1,0,1,1,0,1,0,1;0,1,0,1,1,1,0,1,0;1,0,1,0,1,0,1,0,1;1,0,1,0,1,1,0,1,0];
            trueResForwardFlux = [1,2;1,3;2,4;3,4;4,5;5,6;5,7;6,1;7,1];
            trueResMetList = {'*A','*B','*P','*Q'};
           
            testCase.verifyEqual(trueResRevMatrix, revMatrix);  
            testCase.verifyEqual(trueResForwardFlux, forwardFlux);  
            testCase.verifyEqual(trueResMetList, metList);  
            
        end
    end
end

