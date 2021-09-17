classdef readInputTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testReadInputSubstrateInhibOrderedBiBi(testCase)
            
            filename = 'substrateInhibOrderedBiBi';
                
            [nodeList,edge,kineticMatrix,forwardFlux] = readInput(filename);
            
            trueResNodeList = [1,2,3,4,5,6];
            trueResEdge = {3,[2,3,6],[1,2,3];2,[1,4],[1,4];...
                           1,1,2;2,[2,5],[4,5];2,[4,6],[5,6];...
                           2,[1,5],[3,6]};
            trueResKineticMatrix = {0,'k01.*A','k03.*I',0,0,'k12.*Q';...
                                    'k02',0,0,'k05.*B',0,0;'k04',0,0,0,0,0;...
                                    0,'k06',0,0,'k07',0;0,0,0,'k08',0,'k09';...
                                    'k11',0,0,0,'k10.*P',0};
            trueResForwardFlux = [1,2;1,3;2,4;4,5;5,6;6,1];
            
            testCase.verifyEqual(nodeList, trueResNodeList);
            testCase.verifyEqual(edge, trueResEdge);
            testCase.verifyEqual(kineticMatrix, trueResKineticMatrix);
            testCase.verifyEqual(forwardFlux, trueResForwardFlux);
            
        end
        
        function testReadInputOrdBiBiInhibPromiscuous(testCase)
            
            filename = 'OrdPromiscCompInhibIndep';
                
            [nodeList,edge,kineticMatrix,forwardFlux] = readInput(filename);
            
            trueResNodeList = [1,2,3,4,5,6,7,8,9,10,11];
            trueResEdge = {4,[2,6,7,11],[1,2,3,4];3,[1,3,4],[1,5,6];...
                           1,2,5;2,[2,5],[6,7];2,[4,6],[7,8];2,[1,5],[2,8];...
                           3,[1,8,9],[3,9,10];1,7,9;2,[7,10],[10,11];...
                           2,[9,11],[11,12];2,[1,10],[4,12]};
            trueResKineticMatrix = {0,'k01.*A',0,0,0,'k12.*Q','k13.*C',0,0,0,'k24.*R';...
                                    'k02',0,'k03.*I','k05.*B',0,0,0,0,0,0,0;...
                                    0,'k04',0,0,0,0,0,0,0,0,0;...
                                    0,'k06',0,0,'k07',0,0,0,0,0,0;0,0,0,'k08',0,'k09',0,0,0,0,0;...
                                    'k11',0,0,0,'k10.*P1',0,0,0,0,0,0;...
                                    'k14',0,0,0,0,0,0,'k15.*I','k17.*D',0,0;...
                                    0,0,0,0,0,0,'k16',0,0,0,0;0,0,0,0,0,0,...
                                    'k18',0,0,'k19',0;0,0,0,0,0,0,0,0,'k20',...
                                    0,'k21';'k23',0,0,0,0,0,0,0,0,'k22.*P2',0};
            trueResForwardFlux = [1,2;2,3;2,4;4,5;5,6;6,1;1,7;7,8;7,9;9,10;10,11;11,1];
            
            testCase.verifyEqual(nodeList, trueResNodeList);
            testCase.verifyEqual(edge, trueResEdge);
            testCase.verifyEqual(kineticMatrix, trueResKineticMatrix);
            testCase.verifyEqual(forwardFlux, trueResForwardFlux);
            
        end
        
        function testReadInputRandomBiBi(testCase)
            
            filename = 'randomBiBi';
                
            [nodeList,edge,kineticMatrix,forwardFlux] = readInput(filename);
            
            trueResNodeList = [1,2,3,4,5,6,7];
            trueResEdge = {4,[2,3,6,7],[1,2,3,4];2,[1,4],[1,5];2,[1,4],[2,6];...
                           3,[2,3,5],[5,6,7];3,[4,6,7],[7,8,9];2,[1,5],[3,8];...
                           2,[1,5],[4,9]};
            trueResKineticMatrix = {0,'k01.*A','k03.*B',0,0,'k16.*Q','k18.*P';...
                                    'k02',0,0,'k05.*B',0,0,0;'k04',0,0,'k07.*A',0,0,0;...
                                    0,'k06','k08',0,'k09',0,0;0,0,0,'k10',0,'k11','k13';...
                                    'k15',0,0,0,'k12.*P',0,0;'k17',0,0,0,'k14.*Q',0,0};
            trueResForwardFlux = [1,2;1,3;2,4;3,4;4,5;5,6;5,7;6,1;7,1];
            
            testCase.verifyEqual(nodeList, trueResNodeList);
            testCase.verifyEqual(edge, trueResEdge);
            testCase.verifyEqual(kineticMatrix, trueResKineticMatrix);
            testCase.verifyEqual(forwardFlux, trueResForwardFlux);
            
        end
    end
end

