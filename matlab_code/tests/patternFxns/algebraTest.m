classdef algebraTest < matlab.unittest.TestCase
   
    methods (Test)
        function testAlgebra1(testCase)
                      
            B = [];
            C = [1;4];
            A = algebra(B,C);
            
            trueRes = [1;4];
            
            testCase.verifyEqual(A, trueRes);            
        end
        
        function testAlgebra2(testCase)
                        
            B = [1;4];
            C = [2];
            A = algebra(B,C);
            
            trueRes = [1 2; 2 4];
            
            testCase.verifyEqual(A, trueRes);            
        end
        
        function testAlgebra3(testCase)
                        
            B = [1 2; 2 4];
            C = [4;5];
            A = algebra(B,C);
            
            trueRes = [1 2 4; 1 2 5; 2 4 5];
            
            testCase.verifyEqual(A,trueRes);            
        end
        
        function testAlgebra4(testCase)
                        
            B = [1 2 4; 1 2 5; 2 4 5];
            C = [5;6];
            A = algebra(B,C);
            
            trueRes = [1 2 4 5; 1 2 4 6; 1 2 5 6; 2 4 5 6];
            
            testCase.verifyEqual(A, trueRes);            
        end
    end
end

