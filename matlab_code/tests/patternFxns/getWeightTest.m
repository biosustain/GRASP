classdef getWeightTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testGetWeight1(testCase)
            
            shortPath = [1,2,4,5];
            forwardFlux = [1,2;1,3;2,4;4,5;5,6;6,1];
            
            W = getWeight(shortPath,forwardFlux);
            
            trueResW = [2;1;2;2;1;1];
           
            testCase.verifyEqual(W, trueResW);            
        end

        function testGetWeight2(testCase)
            
            shortPath = [5,6,1,2];
            forwardFlux = [1,2;1,3;2,4;4,5;5,6;6,1];
            
            W = getWeight(shortPath,forwardFlux);
            
            trueResW = [2;1;1;1;2;2];
           
            testCase.verifyEqual(W, trueResW);            
        end

        function testGetWeight3(testCase)
            
            shortPath = [5,1];
            forwardFlux = [1,2;2,3;3,1;1,4;4,5;5,1];
            
            W = getWeight(shortPath,forwardFlux);
            
            trueResW = [1;1;1;1;1;2];
           
            testCase.verifyEqual(W, trueResW);            
        end

        function testGetWeight4(testCase)
            
            shortPath = [10,11,1];
            forwardFlux = [1,2;2,3;2,4;4,5;5,6;6,1;1,7;7,8;7,9;9,10;10,11;11,1];
            
            W = getWeight(shortPath,forwardFlux);
            
            trueResW = [1;1;1;1;1;1;1;1;1;1;2;2];
           
            testCase.verifyEqual(W, trueResW);            
        end

        function testGetWeight5(testCase)
            
            shortPath = [2,4,5];
            forwardFlux = [1,2;1,3;2,4;3,4;4,5;5,6;5,7;6,1;7,1];
            
            W = getWeight(shortPath,forwardFlux);
            
            trueResW = [1;1;2;1;2;1;1;1;1];
           
            testCase.verifyEqual(W, trueResW);            
        end
    end
end

