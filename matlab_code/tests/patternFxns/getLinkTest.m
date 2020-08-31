classdef  getLinkTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testGetLink1(testCase)
            
            edge = {3,[2,3,6],[1,2,3];2,[1,4],[1,4];1,1,2;2,[2,5],[4,5];...
                    2,[4,6],[5,6];2,[1,5],[3,6]};
                
            [edgeMatrix,edgeList] = getLink(edge);
            
            trueResEdgeMatrix = [0,1,2,0,0,3;1,0,0,4,0,0;2,0,0,0,0,0;0,4,0,0,5,0;0,0,0,5,0,6;3,0,0,0,6,0];
            trueResEdgeList = [1,2;1,3;1,6;2,4;4,5;5,6];
           
            testCase.verifyEqual(edgeMatrix, trueResEdgeMatrix);
            testCase.verifyEqual(edgeList, trueResEdgeList);
            
        end
        
        function testGetLinkPromiscuous(testCase)
            
            edge = {4,[2,6,7,11],[1,2,3,4];3,[1,3,4],[1,5,6];1,2,5;2,[2,5],[6,7];...
                    2,[4,6],[7,8];2,[1,5],[2,8];3,[1,8,9],[3,9,10];...
                    1,7,9;2,[7,10],[10,11];2,[9,11],[11,12];2,[1,10],[4,12]};
                
            [edgeMatrix,edgeList] = getLink(edge);
            
            trueResEdgeMatrix = [0,1,0,0,0,2,3,0,0,0,4;1,0,5,6,0,0,0,0,0,0,0;...
                                 0,5,0,0,0,0,0,0,0,0,0;0,6,0,0,7,0,0,0,0,0,0;...
                                 0,0,0,7,0,8,0,0,0,0,0;2,0,0,0,8,0,0,0,0,0,0;...
                                 3,0,0,0,0,0,0,9,10,0,0;0,0,0,0,0,0,9,0,0,0,0;...
                                 0,0,0,0,0,0,10,0,0,11,0;0,0,0,0,0,0,0,0,11,0,12;...
                                 4,0,0,0,0,0,0,0,0,12,0];
            trueResEdgeList = [1,2;1,6;1,7;1,11;2,3;2,4;4,5;5,6;7,8;7,9;9,10;10,11];
           
            testCase.verifyEqual(edgeMatrix, trueResEdgeMatrix);
            testCase.verifyEqual(edgeList, trueResEdgeList);
            
        end
        
        function testGetLinkRandomBiBi(testCase)
            
            edge = {4,[2,3,6,7],[1,2,3,4];2,[1,4],[1,5];2,[1,4],[2,6];...
                    3,[2,3,5],[5,6,7];3,[4,6,7],[7,8,9];2,[1,5],[3,8];...
                    2,[1,5],[4,9]};
                
            [edgeMatrix,edgeList] = getLink(edge);
            
            trueResEdgeMatrix = [0,1,2,0,0,3,4;1,0,0,5,0,0,0;2,0,0,6,0,0,0;...
                                 0,5,6,0,7,0,0;0,0,0,7,0,8,9;3,0,0,0,8,0,0;...
                                 4,0,0,0,9,0,0];
            trueResEdgeList = [1,2;1,3;1,6;1,7;2,4;3,4;4,5;5,6;5,7];
           
            testCase.verifyEqual(edgeMatrix, trueResEdgeMatrix);
            testCase.verifyEqual(edgeList, trueResEdgeList);
            
        end
    end
end

