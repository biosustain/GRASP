classdef colVectorTest < matlab.unittest.TestCase
    
    methods (Test)
        function testColVector1(testCase)
            
            vector = {[1,2,4,5,6,1];[2,4,5,6,1,2]};            
            vector = colVector(vector);
            
            trueRes = {[1,2,4,5,6,1];[2,4,5,6,1,2]};
           
            testCase.verifyEqual(trueRes,vector);
            
        end
        
        function testColVector2(testCase)
            
            vector = {[1,2,3,1],[1,4,5,1]};           
            vector = colVector(vector);
            
            trueRes = {[1,2,3,1];[1,4,5,1]};
           
            testCase.verifyEqual(trueRes,vector);
            
        end
    end
end

