classdef cycledPathsTest < matlab.unittest.TestCase
    
    methods (Test)
        function testCycledPaths1(testCase)
            subsNodes = [1;2];
            prodNodes = [5;6];
            forwardFlux = [1,2;1,3;2,4;4,5;5,6;6,1];
            paths = cycledPaths(subsNodes,prodNodes,forwardFlux);
            
            trueRes = {[1,2,4,5,6,1];[2,4,5,6,1,2]};
           
            testCase.verifyEqual(paths, trueRes);
            
        end
        
        function testCycledPaths2(testCase)
            subsNodes = [1];
            prodNodes = [3;5];
            forwardFlux = [1,2;2,3;3,1;1,4;4,5;5,1];
            paths = cycledPaths(subsNodes,prodNodes,forwardFlux);
            
            trueRes = {[1,2,3,1],[1,4,5,1]};
           
            testCase.verifyEqual(paths, trueRes);
            
        end
        
        function testCycledPaths3(testCase)
            subsNodes = [1;2;7];
            prodNodes = [5;6;10;11];
            forwardFlux = [1,2;2,3;2,4;4,5;5,6;6,1;1,7;7,8;7,9;9,10;10,11;11,1];
            paths = cycledPaths(subsNodes,prodNodes,forwardFlux);
            
            trueRes = {[1,2,4,5,6,1],[1,7,9,10,11,1];[2,4,5,6,1,2],[2,4,5,6,1,2];[7,9,10,11,1,7],[7,9,10,11,1,7]};
           
            testCase.verifyEqual(paths, trueRes);
            
        end
        
        function testCycledPaths4(testCase)
            subsNodes = [1;2;3];
            prodNodes = [5;6;7];
            forwardFlux = [1,2;1,3;2,4;3,4;4,5;5,6;5,7;6,1;7,1];
            paths = cycledPaths(subsNodes,prodNodes,forwardFlux);
            
            trueRes = {[1,2,4,5,6,1];[2,4,5,6,1,2];[3,4,5,6,1,3];[1,3,4,5,6,1];...
                       [1,3,4,5,7,1];[1,2,4,5,7,1];[2,4,5,7,1,2];[3,4,5,7,1,3]};
           
            testCase.verifyEqual(paths, trueRes);
            
        end
     
    end
end

