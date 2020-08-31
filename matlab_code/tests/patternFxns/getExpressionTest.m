classdef getExpressionTest < matlab.unittest.TestCase
    
    methods (Test)
        function testGetExpression1(testCase)
            
            Nodes = [1,2,3,4,5,6];
            PatternNumber = 5;
        	pathway(:, :, 1) = [1, 2;...
                                1, 3;...
                                1, 6;...
                                2, 4;...
                                4, 5];
            pathway(:, :, 2) = [1, 2;...
                                1, 3;...
                                1, 6;...
                                2, 4;...
                                5, 6];
            pathway(:, :, 3) = [1, 2;...
                                1, 3;...
                                1, 6;...
                                4, 5;...
                                5, 6];
            pathway(:, :, 4) = [1, 2;...
                                1, 3;...
                                2, 4;...
                                4, 5;...
                                5, 6];
            pathway(:, :, 5) = [1, 3;...
                                1, 6;...
                                2, 4;...
                                4, 5;...
                                5, 6];                
            KineticMatrix = {0,'k01.*A','k03.*I',0,0,'k12.*Q';'k02',0,0,'k05.*B',0,0;...
                            'k04',0,0,0,0,0;0,'k06',0,0,'k07',0;0,0,0,'k08',0,'k09';...
                            'k11',0,0,0,'k10.*P',0};
            expression = cell(length(Nodes),PatternNumber);
            
            for i = 1:length(Nodes)
                node = Nodes(i);
                for j = 1:PatternNumber
                    exp = '';
                    exp = getExpression(node,pathway(:,:,j),KineticMatrix,exp);
                    expression{i,j} = exp;
                end
            end
            
            trueRes = {'k02*k04*k11*k06*k08','k02*k04*k11*k06*k09','k02*k04*k11*k09*k07','k02*k04*k06*k08*k10.*P','k04*k11*k09*k07*k05.*B';...
                       'k01.*A*k06*k04*k11*k08','k01.*A*k06*k04*k11*k09','k01.*A*k04*k11*k09*k07','k01.*A*k06*k04*k08*k10.*P','k06*k08*k10.*P*k12.*Q*k04';...
                       'k03.*I*k02*k11*k06*k08','k03.*I*k02*k11*k06*k09','k03.*I*k02*k11*k09*k07','k03.*I*k02*k06*k08*k10.*P','k03.*I*k11*k09*k07*k05.*B';...
                       'k05.*B*k08*k01.*A*k04*k11','k05.*B*k01.*A*k04*k11*k09','k08*k10.*P*k12.*Q*k02*k04','k05.*B*k08*k01.*A*k10.*P*k04','k05.*B*k08*k10.*P*k12.*Q*k04';...
                       'k07*k05.*B*k01.*A*k04*k11','k10.*P*k12.*Q*k02*k04*k06','k07*k10.*P*k12.*Q*k02*k04','k07*k10.*P*k05.*B*k01.*A*k04','k07*k10.*P*k05.*B*k12.*Q*k04';...
                       'k12.*Q*k02*k04*k06*k08','k12.*Q*k09*k02*k04*k06','k12.*Q*k09*k02*k04*k07','k09*k07*k05.*B*k01.*A*k04','k12.*Q*k09*k04*k07*k05.*B'};
           
            testCase.verifyEqual(expression, trueRes);
            
        end
        
        function testGetExpressionPromiscuous(testCase)
            
            Nodes = [1,2,3,4,5,6,7,8,9,10,11];
            PatternNumber = 25;
            currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');                      
        	pathway = load(fullfile(currentPath{1}, 'testFiles', 'pathwayPromiscuous.mat'));   
            pathway = pathway.pathway;
            KineticMatrix = {0,'k01.*A',0,0,0,'k12.*Q','k13.*C',0,0,0,'k24.*R';...
                             'k02',0,'k03.*I','k05.*B',0,0,0,0,0,0,0;0,'k04',0,0,0,0,0,0,0,0,0;...
                             0,'k06',0,0,'k07',0,0,0,0,0,0;0,0,0,'k08',0,'k09',0,0,0,0,0;...
                             'k11',0,0,0,'k10.*P1',0,0,0,0,0,0;'k14',0,0,0,0,0,0,'k15.*I','k17.*D',0,0;...
                             0,0,0,0,0,0,'k16',0,0,0,0;0,0,0,0,0,0,'k18',0,0,'k19',0;...
                             0,0,0,0,0,0,0,0,'k20',0,'k21';'k23',0,0,0,0,0,0,0,0,'k22.*P2',0};
            expression = cell(length(Nodes),PatternNumber);
            
            for i = 1:length(Nodes)
                node = Nodes(i);
                for j = 1:PatternNumber
                    exp = '';
                    exp = getExpression(node,pathway(:,:,j),KineticMatrix,exp);
                    expression{i,j} = exp;
                end
            end
           
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResExpressionPromiscuous.mat'));
            trueRes = trueRes.trueRes;
            
            testCase.verifyEqual(expression, trueRes);
            
        end
        
        function testGetExpressionRandomBiBi(testCase)
            
            Nodes = [1,2,3,4,5,6,7];
            PatternNumber = 48;
            currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');                      
        	pathway = load(fullfile(currentPath{1}, 'testFiles', 'pathwayRandomBiBi.mat'));   
            pathway = pathway.pathway;
            KineticMatrix = {0,'k01.*A','k03.*B',0,0,'k16.*Q','k18.*P';'k02',0,0,'k05.*B',0,0,0;...
                             'k04',0,0,'k07.*A',0,0,0;0,'k06','k08',0,'k09',0,0;0,0,0,'k10',0,'k11','k13';...
                             'k15',0,0,0,'k12.*P',0,0;'k17',0,0,0,'k14.*Q',0,0};
            expression = cell(length(Nodes),PatternNumber);
            
            for i = 1:length(Nodes)
                node = Nodes(i);
                for j = 1:PatternNumber
                    exp = '';
                    exp = getExpression(node,pathway(:,:,j),KineticMatrix,exp);
                    expression{i,j} = exp;
                end
            end
            
            trueRes = load(fullfile(currentPath{1}, 'testFiles', 'trueResExpressionRandomBiBi.mat'));    
            trueRes = trueRes.expression;
           
            testCase.verifyEqual(expression, trueRes);
            
        end
        
    end
end

