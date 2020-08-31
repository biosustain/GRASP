classdef buildReactionTest < matlab.unittest.TestCase
    
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
        function testBuildReaction(testCase)
            state = {'E1= k02*k04*k11*k06*k08+k02*k04*k11*k06*k09+k02*k04*k11*k09*k07+k02*k04*k06*k08*k10.*P+k04*k11*k09*k07*k05.*B;', ...
                     'E2= k01.*A*k06*k04*k11*k08+k01.*A*k06*k04*k11*k09+k01.*A*k04*k11*k09*k07+k01.*A*k06*k04*k08*k10.*P+k06*k08*k10.*P*k12.*Q*k04;', ...
                     'E3= k03.*I*k02*k11*k06*k08+k03.*I*k02*k11*k06*k09+k03.*I*k02*k11*k09*k07+k03.*I*k02*k06*k08*k10.*P+k03.*I*k11*k09*k07*k05.*B;', ...
                     'E4= k05.*B*k08*k01.*A*k04*k11+k05.*B*k01.*A*k04*k11*k09+k08*k10.*P*k12.*Q*k02*k04+k05.*B*k08*k01.*A*k10.*P*k04+k05.*B*k08*k10.*P*k12.*Q*k04;', ...
                     'E5= k07*k05.*B*k01.*A*k04*k11+k10.*P*k12.*Q*k02*k04*k06+k07*k10.*P*k12.*Q*k02*k04+k07*k10.*P*k05.*B*k01.*A*k04+k07*k10.*P*k05.*B*k12.*Q*k04;', ...
                     'E6= k12.*Q*k02*k04*k06*k08+k12.*Q*k09*k02*k04*k06+k12.*Q*k09*k02*k04*k07+k09*k07*k05.*B*k01.*A*k04+k12.*Q*k09*k04*k07*k05.*B;'};
            rateList = {'k01', 'k02', 'k03', 'k04', 'k05', 'k06', 'k07', 'k08', ...
                        'k09', 'k10', 'k11', 'k12'};
            metList = {'*A', '*B', '*I', '*P', '*Q'};
            numTerm = {'+k09', '-k10.*P'};
            prodNum = [5 6];
            reactionName = 'r_r11';
            promiscuousRxnI = 0;
            
            buildReaction(state,rateList,metList,numTerm,prodNum,reactionName, promiscuousRxnI)
           
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildReaction.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);   
         end
        
        function testBuildReactionPromiscuous1(testCase)

           state = {'E1= k02*k11*k14*k23*k04*k06*k16*k18*k08*k20+k02*k11*k14*k23*k04*k06*k16*k18*k21*k08+k02*k11*k14*k23*k04*k06*k16*k21*k08*k19+k02*k11*k14*k23*k04*k06*k09*k16*k18*k20+k02*k11*k14*k23*k04*k06*k09*k16*k18*k21+k02*k11*k14*k23*k04*k06*k09*k16*k21*k19+k02*k11*k14*k23*k04*k09*k16*k18*k07*k20+k02*k11*k14*k23*k04*k09*k16*k18*k21*k07+k02*k11*k14*k23*k04*k09*k16*k21*k07*k19+k02*k11*k14*k04*k06*k16*k18*k08*k20*k22.*P2+k02*k11*k14*k04*k06*k09*k16*k18*k20*k22.*P2+k02*k11*k14*k04*k09*k16*k18*k07*k20*k22.*P2+k02*k11*k23*k04*k06*k21*k08*k19*k17.*D*k16+k02*k11*k23*k04*k06*k09*k21*k19*k17.*D*k16+k02*k11*k23*k04*k09*k21*k07*k19*k17.*D*k16+k02*k14*k23*k04*k06*k16*k18*k08*k20*k10.*P1+k02*k14*k23*k04*k06*k16*k18*k21*k08*k10.*P1+k02*k14*k23*k04*k06*k16*k21*k08*k19*k10.*P1+k02*k14*k04*k06*k16*k18*k08*k20*k10.*P1*k22.*P2+k02*k23*k04*k06*k21*k08*k19*k10.*P1*k17.*D*k16+k11*k14*k23*k09*k16*k18*k07*k20*k05.*B*k04+k11*k14*k23*k09*k16*k18*k21*k07*k05.*B*k04+k11*k14*k23*k09*k16*k21*k07*k19*k05.*B*k04+k11*k14*k09*k16*k18*k07*k20*k05.*B*k22.*P2*k04+k11*k23*k09*k21*k07*k19*k05.*B*k17.*D*k04*k16;', ...
                    'E2= k01.*A*k04*k06*k11*k14*k23*k08*k16*k18*k20+k01.*A*k04*k06*k11*k14*k23*k08*k16*k18*k21+k01.*A*k04*k06*k11*k14*k23*k08*k16*k21*k19+k01.*A*k04*k06*k11*k14*k23*k09*k16*k18*k20+k01.*A*k04*k06*k11*k14*k23*k09*k16*k18*k21+k01.*A*k04*k06*k11*k14*k23*k09*k16*k21*k19+k01.*A*k04*k11*k14*k23*k09*k16*k18*k07*k20+k01.*A*k04*k11*k14*k23*k09*k16*k18*k21*k07+k01.*A*k04*k11*k14*k23*k09*k16*k21*k07*k19+k01.*A*k04*k06*k11*k14*k08*k16*k18*k20*k22.*P2+k01.*A*k04*k06*k11*k14*k09*k16*k18*k20*k22.*P2+k01.*A*k04*k11*k14*k09*k16*k18*k07*k20*k22.*P2+k01.*A*k04*k06*k11*k23*k08*k21*k19*k17.*D*k16+k01.*A*k04*k06*k11*k23*k09*k21*k19*k17.*D*k16+k01.*A*k04*k11*k23*k09*k21*k07*k19*k17.*D*k16+k01.*A*k04*k06*k14*k23*k08*k16*k18*k10.*P1*k20+k01.*A*k04*k06*k14*k23*k08*k16*k18*k21*k10.*P1+k01.*A*k04*k06*k14*k23*k08*k16*k21*k10.*P1*k19+k01.*A*k04*k06*k14*k08*k16*k18*k10.*P1*k20*k22.*P2+k01.*A*k04*k06*k23*k08*k21*k10.*P1*k19*k17.*D*k16+k04*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k18*k20+k04*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k18*k21+k04*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k21*k19+k04*k06*k08*k10.*P1*k12.*Q*k14*k16*k18*k20*k22.*P2+k04*k06*k08*k10.*P1*k12.*Q*k23*k21*k19*k17.*D*k16;', ...
                    'E3= k03.*I*k01.*A*k06*k11*k14*k23*k08*k16*k18*k20+k03.*I*k01.*A*k06*k11*k14*k23*k08*k16*k18*k21+k03.*I*k01.*A*k06*k11*k14*k23*k08*k16*k21*k19+k03.*I*k01.*A*k06*k11*k14*k23*k09*k16*k18*k20+k03.*I*k01.*A*k06*k11*k14*k23*k09*k16*k18*k21+k03.*I*k01.*A*k06*k11*k14*k23*k09*k16*k21*k19+k03.*I*k01.*A*k11*k14*k23*k09*k16*k18*k07*k20+k03.*I*k01.*A*k11*k14*k23*k09*k16*k18*k21*k07+k03.*I*k01.*A*k11*k14*k23*k09*k16*k21*k07*k19+k03.*I*k01.*A*k06*k11*k14*k08*k16*k18*k20*k22.*P2+k03.*I*k01.*A*k06*k11*k14*k09*k16*k18*k20*k22.*P2+k03.*I*k01.*A*k11*k14*k09*k16*k18*k07*k20*k22.*P2+k03.*I*k01.*A*k06*k11*k23*k08*k21*k19*k17.*D*k16+k03.*I*k01.*A*k06*k11*k23*k09*k21*k19*k17.*D*k16+k03.*I*k01.*A*k11*k23*k09*k21*k07*k19*k17.*D*k16+k03.*I*k01.*A*k06*k14*k23*k08*k16*k18*k10.*P1*k20+k03.*I*k01.*A*k06*k14*k23*k08*k16*k18*k21*k10.*P1+k03.*I*k01.*A*k06*k14*k23*k08*k16*k21*k10.*P1*k19+k03.*I*k01.*A*k06*k14*k08*k16*k18*k10.*P1*k20*k22.*P2+k03.*I*k01.*A*k06*k23*k08*k21*k10.*P1*k19*k17.*D*k16+k03.*I*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k18*k20+k03.*I*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k18*k21+k03.*I*k06*k08*k10.*P1*k12.*Q*k14*k23*k16*k21*k19+k03.*I*k06*k08*k10.*P1*k12.*Q*k14*k16*k18*k20*k22.*P2+k03.*I*k06*k08*k10.*P1*k12.*Q*k23*k21*k19*k17.*D*k16;', ...
                    'E4= k05.*B*k08*k01.*A*k04*k11*k14*k23*k16*k18*k20+k05.*B*k08*k01.*A*k04*k11*k14*k23*k16*k18*k21+k05.*B*k08*k01.*A*k04*k11*k14*k23*k16*k21*k19+k05.*B*k01.*A*k04*k11*k14*k23*k09*k16*k18*k20+k05.*B*k01.*A*k04*k11*k14*k23*k09*k16*k18*k21+k05.*B*k01.*A*k04*k11*k14*k23*k09*k16*k21*k19+k08*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k18*k20+k08*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k18*k21+k08*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k21*k19+k05.*B*k08*k01.*A*k04*k11*k14*k16*k18*k20*k22.*P2+k05.*B*k01.*A*k04*k11*k14*k09*k16*k18*k20*k22.*P2+k08*k10.*P1*k12.*Q*k02*k14*k04*k16*k18*k20*k22.*P2+k05.*B*k08*k01.*A*k04*k11*k23*k21*k19*k17.*D*k16+k05.*B*k01.*A*k04*k11*k23*k09*k21*k19*k17.*D*k16+k08*k10.*P1*k12.*Q*k02*k23*k04*k21*k19*k17.*D*k16+k05.*B*k08*k01.*A*k04*k10.*P1*k14*k23*k16*k18*k20+k05.*B*k08*k01.*A*k04*k10.*P1*k14*k23*k16*k18*k21+k05.*B*k08*k01.*A*k04*k10.*P1*k14*k23*k16*k21*k19+k05.*B*k08*k01.*A*k04*k10.*P1*k14*k16*k18*k20*k22.*P2+k05.*B*k08*k01.*A*k04*k10.*P1*k23*k21*k19*k17.*D*k16+k05.*B*k08*k04*k10.*P1*k12.*Q*k14*k23*k16*k18*k20+k05.*B*k08*k04*k10.*P1*k12.*Q*k14*k23*k16*k18*k21+k05.*B*k08*k04*k10.*P1*k12.*Q*k14*k23*k16*k21*k19+k05.*B*k08*k04*k10.*P1*k12.*Q*k14*k16*k18*k20*k22.*P2+k05.*B*k08*k04*k10.*P1*k12.*Q*k23*k21*k19*k17.*D*k16;', ...
                    'E5= k07*k05.*B*k01.*A*k04*k11*k14*k23*k16*k18*k20+k07*k05.*B*k01.*A*k04*k11*k14*k23*k16*k18*k21+k07*k05.*B*k01.*A*k04*k11*k14*k23*k16*k21*k19+k10.*P1*k12.*Q*k02*k14*k23*k04*k06*k16*k18*k20+k10.*P1*k12.*Q*k02*k14*k23*k04*k06*k16*k18*k21+k10.*P1*k12.*Q*k02*k14*k23*k04*k06*k16*k21*k19+k07*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k18*k20+k07*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k18*k21+k07*k10.*P1*k12.*Q*k02*k14*k23*k04*k16*k21*k19+k07*k05.*B*k01.*A*k04*k11*k14*k16*k18*k20*k22.*P2+k10.*P1*k12.*Q*k02*k14*k04*k06*k16*k18*k20*k22.*P2+k07*k10.*P1*k12.*Q*k02*k14*k04*k16*k18*k20*k22.*P2+k07*k05.*B*k01.*A*k04*k11*k23*k21*k19*k17.*D*k16+k10.*P1*k12.*Q*k02*k23*k04*k06*k21*k19*k17.*D*k16+k07*k10.*P1*k12.*Q*k02*k23*k04*k21*k19*k17.*D*k16+k07*k10.*P1*k05.*B*k01.*A*k04*k14*k23*k16*k18*k20+k07*k10.*P1*k05.*B*k01.*A*k04*k14*k23*k16*k18*k21+k07*k10.*P1*k05.*B*k01.*A*k04*k14*k23*k16*k21*k19+k07*k10.*P1*k05.*B*k01.*A*k04*k14*k16*k18*k20*k22.*P2+k07*k10.*P1*k05.*B*k01.*A*k04*k23*k21*k19*k17.*D*k16+k07*k10.*P1*k05.*B*k12.*Q*k04*k14*k23*k16*k18*k20+k07*k10.*P1*k05.*B*k12.*Q*k04*k14*k23*k16*k18*k21+k07*k10.*P1*k05.*B*k12.*Q*k04*k14*k23*k16*k21*k19+k07*k10.*P1*k05.*B*k12.*Q*k04*k14*k16*k18*k20*k22.*P2+k07*k10.*P1*k05.*B*k12.*Q*k04*k23*k21*k19*k17.*D*k16;', ...
                    'E6= k12.*Q*k02*k14*k23*k04*k06*k16*k18*k08*k20+k12.*Q*k02*k14*k23*k04*k06*k16*k18*k21*k08+k12.*Q*k02*k14*k23*k04*k06*k16*k21*k08*k19+k12.*Q*k09*k02*k14*k23*k04*k06*k16*k18*k20+k12.*Q*k09*k02*k14*k23*k04*k06*k16*k18*k21+k12.*Q*k09*k02*k14*k23*k04*k06*k16*k21*k19+k12.*Q*k09*k02*k14*k23*k07*k04*k16*k18*k20+k12.*Q*k09*k02*k14*k23*k07*k04*k16*k18*k21+k12.*Q*k09*k02*k14*k23*k07*k04*k16*k21*k19+k12.*Q*k02*k14*k04*k06*k16*k18*k08*k20*k22.*P2+k12.*Q*k09*k02*k14*k04*k06*k16*k18*k20*k22.*P2+k12.*Q*k09*k02*k14*k07*k04*k16*k18*k20*k22.*P2+k12.*Q*k02*k23*k04*k06*k21*k08*k19*k17.*D*k16+k12.*Q*k09*k02*k23*k04*k06*k21*k19*k17.*D*k16+k12.*Q*k09*k02*k23*k07*k04*k21*k19*k17.*D*k16+k09*k07*k05.*B*k01.*A*k04*k14*k23*k16*k18*k20+k09*k07*k05.*B*k01.*A*k04*k14*k23*k16*k18*k21+k09*k07*k05.*B*k01.*A*k04*k14*k23*k16*k21*k19+k09*k07*k05.*B*k01.*A*k04*k14*k16*k18*k20*k22.*P2+k09*k07*k05.*B*k01.*A*k04*k23*k21*k19*k17.*D*k16+k12.*Q*k09*k14*k23*k07*k16*k18*k05.*B*k20*k04+k12.*Q*k09*k14*k23*k07*k16*k18*k21*k05.*B*k04+k12.*Q*k09*k14*k23*k07*k16*k21*k05.*B*k19*k04+k12.*Q*k09*k14*k07*k16*k18*k05.*B*k20*k04*k22.*P2+k12.*Q*k09*k23*k07*k21*k05.*B*k19*k04*k17.*D*k16;', ...
                    'E7= k13.*C*k16*k18*k02*k11*k23*k20*k04*k06*k08+k13.*C*k16*k18*k02*k11*k23*k04*k06*k21*k08+k13.*C*k16*k02*k11*k23*k04*k06*k21*k08*k19+k13.*C*k16*k18*k02*k11*k23*k20*k04*k06*k09+k13.*C*k16*k18*k02*k11*k23*k04*k06*k09*k21+k13.*C*k16*k02*k11*k23*k04*k06*k09*k21*k19+k13.*C*k16*k18*k02*k11*k23*k20*k04*k09*k07+k13.*C*k16*k18*k02*k11*k23*k04*k09*k21*k07+k13.*C*k16*k02*k11*k23*k04*k09*k21*k07*k19+k13.*C*k16*k18*k02*k11*k20*k04*k06*k22.*P2*k08+k13.*C*k16*k18*k02*k11*k20*k04*k06*k09*k22.*P2+k13.*C*k16*k18*k02*k11*k20*k04*k09*k22.*P2*k07+k16*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k06*k08+k16*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k06*k09+k16*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k09*k07+k13.*C*k16*k18*k02*k23*k20*k04*k06*k08*k10.*P1+k13.*C*k16*k18*k02*k23*k04*k06*k21*k08*k10.*P1+k13.*C*k16*k02*k23*k04*k06*k21*k08*k19*k10.*P1+k13.*C*k16*k18*k02*k20*k04*k06*k22.*P2*k08*k10.*P1+k16*k18*k20*k22.*P2*k24.*R*k02*k04*k06*k08*k10.*P1+k13.*C*k16*k18*k11*k23*k20*k09*k07*k05.*B*k04+k13.*C*k16*k18*k11*k23*k09*k21*k07*k05.*B*k04+k13.*C*k16*k11*k23*k09*k21*k07*k19*k05.*B*k04+k13.*C*k16*k18*k11*k20*k09*k22.*P2*k07*k05.*B*k04+k16*k18*k20*k22.*P2*k24.*R*k11*k09*k07*k05.*B*k04;', ...
                    'E8= k15.*I*k13.*C*k18*k02*k11*k23*k20*k04*k06*k08+k15.*I*k13.*C*k18*k02*k11*k23*k04*k06*k21*k08+k15.*I*k13.*C*k02*k11*k23*k04*k06*k21*k08*k19+k15.*I*k13.*C*k18*k02*k11*k23*k20*k04*k06*k09+k15.*I*k13.*C*k18*k02*k11*k23*k04*k06*k09*k21+k15.*I*k13.*C*k02*k11*k23*k04*k06*k09*k21*k19+k15.*I*k13.*C*k18*k02*k11*k23*k20*k04*k09*k07+k15.*I*k13.*C*k18*k02*k11*k23*k04*k09*k21*k07+k15.*I*k13.*C*k02*k11*k23*k04*k09*k21*k07*k19+k15.*I*k13.*C*k18*k02*k11*k20*k04*k06*k22.*P2*k08+k15.*I*k13.*C*k18*k02*k11*k20*k04*k06*k09*k22.*P2+k15.*I*k13.*C*k18*k02*k11*k20*k04*k09*k22.*P2*k07+k15.*I*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k06*k08+k15.*I*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k06*k09+k15.*I*k18*k20*k22.*P2*k24.*R*k02*k11*k04*k09*k07+k15.*I*k13.*C*k18*k02*k23*k20*k04*k06*k08*k10.*P1+k15.*I*k13.*C*k18*k02*k23*k04*k06*k21*k08*k10.*P1+k15.*I*k13.*C*k02*k23*k04*k06*k21*k08*k19*k10.*P1+k15.*I*k13.*C*k18*k02*k20*k04*k06*k22.*P2*k08*k10.*P1+k15.*I*k18*k20*k22.*P2*k24.*R*k02*k04*k06*k08*k10.*P1+k15.*I*k13.*C*k18*k11*k23*k20*k09*k07*k05.*B*k04+k15.*I*k13.*C*k18*k11*k23*k09*k21*k07*k05.*B*k04+k15.*I*k13.*C*k11*k23*k09*k21*k07*k19*k05.*B*k04+k15.*I*k13.*C*k18*k11*k20*k09*k22.*P2*k07*k05.*B*k04+k15.*I*k18*k20*k22.*P2*k24.*R*k11*k09*k07*k05.*B*k04;', ...
                    'E9= k17.*D*k20*k13.*C*k16*k02*k11*k23*k04*k06*k08+k17.*D*k13.*C*k16*k02*k11*k23*k04*k06*k21*k08+k20*k22.*P2*k24.*R*k02*k11*k14*k04*k06*k16*k08+k17.*D*k20*k13.*C*k16*k02*k11*k23*k04*k06*k09+k17.*D*k13.*C*k16*k02*k11*k23*k04*k06*k09*k21+k20*k22.*P2*k24.*R*k02*k11*k14*k04*k06*k09*k16+k17.*D*k20*k13.*C*k16*k02*k11*k23*k04*k09*k07+k17.*D*k13.*C*k16*k02*k11*k23*k04*k09*k21*k07+k20*k22.*P2*k24.*R*k02*k11*k14*k04*k09*k16*k07+k17.*D*k20*k13.*C*k16*k22.*P2*k02*k11*k04*k06*k08+k17.*D*k20*k13.*C*k16*k22.*P2*k02*k11*k04*k06*k09+k17.*D*k20*k13.*C*k16*k22.*P2*k02*k11*k04*k09*k07+k17.*D*k20*k16*k22.*P2*k24.*R*k02*k11*k04*k06*k08+k17.*D*k20*k16*k22.*P2*k24.*R*k02*k11*k04*k06*k09+k17.*D*k20*k16*k22.*P2*k24.*R*k02*k11*k04*k09*k07+k17.*D*k20*k13.*C*k16*k02*k23*k04*k06*k08*k10.*P1+k17.*D*k13.*C*k16*k02*k23*k04*k06*k21*k08*k10.*P1+k20*k22.*P2*k24.*R*k02*k14*k04*k06*k16*k08*k10.*P1+k17.*D*k20*k13.*C*k16*k22.*P2*k02*k04*k06*k08*k10.*P1+k17.*D*k20*k16*k22.*P2*k24.*R*k02*k04*k06*k08*k10.*P1+k17.*D*k20*k13.*C*k16*k11*k23*k09*k07*k05.*B*k04+k17.*D*k13.*C*k16*k11*k23*k09*k21*k07*k05.*B*k04+k20*k22.*P2*k24.*R*k11*k14*k09*k16*k07*k05.*B*k04+k17.*D*k20*k13.*C*k16*k22.*P2*k11*k09*k07*k05.*B*k04+k17.*D*k20*k16*k22.*P2*k24.*R*k11*k09*k07*k05.*B*k04;', ...
                    'E10= k19*k17.*D*k13.*C*k16*k02*k11*k23*k04*k06*k08+k22.*P2*k24.*R*k02*k11*k14*k04*k06*k16*k18*k08+k19*k22.*P2*k24.*R*k02*k11*k14*k04*k06*k16*k08+k19*k17.*D*k13.*C*k16*k02*k11*k23*k04*k06*k09+k22.*P2*k24.*R*k02*k11*k14*k04*k06*k09*k16*k18+k19*k22.*P2*k24.*R*k02*k11*k14*k04*k06*k09*k16+k19*k17.*D*k13.*C*k16*k02*k11*k23*k04*k09*k07+k22.*P2*k24.*R*k02*k11*k14*k04*k09*k16*k18*k07+k19*k22.*P2*k24.*R*k02*k11*k14*k04*k09*k16*k07+k19*k22.*P2*k17.*D*k13.*C*k16*k02*k11*k04*k06*k08+k19*k22.*P2*k17.*D*k13.*C*k16*k02*k11*k04*k06*k09+k19*k22.*P2*k17.*D*k13.*C*k16*k02*k11*k04*k09*k07+k19*k22.*P2*k17.*D*k24.*R*k16*k02*k11*k04*k06*k08+k19*k22.*P2*k17.*D*k24.*R*k16*k02*k11*k04*k06*k09+k19*k22.*P2*k17.*D*k24.*R*k16*k02*k11*k04*k09*k07+k19*k17.*D*k13.*C*k16*k02*k23*k04*k06*k08*k10.*P1+k22.*P2*k24.*R*k02*k14*k04*k06*k16*k18*k08*k10.*P1+k19*k22.*P2*k24.*R*k02*k14*k04*k06*k16*k08*k10.*P1+k19*k22.*P2*k17.*D*k13.*C*k16*k02*k04*k06*k08*k10.*P1+k19*k22.*P2*k17.*D*k24.*R*k16*k02*k04*k06*k08*k10.*P1+k19*k17.*D*k13.*C*k16*k11*k23*k09*k07*k05.*B*k04+k22.*P2*k24.*R*k11*k14*k09*k16*k18*k07*k05.*B*k04+k19*k22.*P2*k24.*R*k11*k14*k09*k16*k07*k05.*B*k04+k19*k22.*P2*k17.*D*k13.*C*k16*k11*k09*k07*k05.*B*k04+k19*k22.*P2*k17.*D*k24.*R*k16*k11*k09*k07*k05.*B*k04;', ...
                    'E11= k24.*R*k02*k11*k14*k04*k06*k16*k18*k08*k20+k24.*R*k21*k02*k11*k14*k04*k06*k16*k18*k08+k24.*R*k21*k02*k11*k14*k19*k04*k06*k16*k08+k24.*R*k02*k11*k14*k04*k06*k09*k16*k18*k20+k24.*R*k21*k02*k11*k14*k04*k06*k09*k16*k18+k24.*R*k21*k02*k11*k14*k19*k04*k06*k09*k16+k24.*R*k02*k11*k14*k04*k09*k16*k18*k07*k20+k24.*R*k21*k02*k11*k14*k04*k09*k16*k18*k07+k24.*R*k21*k02*k11*k14*k19*k04*k09*k16*k07+k21*k19*k17.*D*k13.*C*k16*k02*k11*k04*k06*k08+k21*k19*k17.*D*k13.*C*k16*k02*k11*k04*k06*k09+k21*k19*k17.*D*k13.*C*k16*k02*k11*k04*k09*k07+k24.*R*k21*k02*k11*k19*k04*k06*k17.*D*k08*k16+k24.*R*k21*k02*k11*k19*k04*k06*k09*k17.*D*k16+k24.*R*k21*k02*k11*k19*k04*k09*k17.*D*k07*k16+k24.*R*k02*k14*k04*k06*k16*k18*k08*k20*k10.*P1+k24.*R*k21*k02*k14*k04*k06*k16*k18*k08*k10.*P1+k24.*R*k21*k02*k14*k19*k04*k06*k16*k08*k10.*P1+k21*k19*k17.*D*k13.*C*k16*k02*k04*k06*k08*k10.*P1+k24.*R*k21*k02*k19*k04*k06*k17.*D*k08*k16*k10.*P1+k24.*R*k11*k14*k09*k16*k18*k07*k20*k05.*B*k04+k24.*R*k21*k11*k14*k09*k16*k18*k07*k05.*B*k04+k24.*R*k21*k11*k14*k19*k09*k16*k07*k05.*B*k04+k21*k19*k17.*D*k13.*C*k16*k11*k09*k07*k05.*B*k04+k24.*R*k21*k11*k19*k09*k17.*D*k07*k16*k05.*B*k04;'};

            rateList = {'k01', 'k02', 'k03', 'k04', 'k05', 'k06', 'k07', 'k08', ...
                        'k09', 'k10', 'k11', 'k12', 'k13', 'k14', 'k15', 'k16',...
                        'k17', 'k18', 'k19', 'k20', 'k21', 'k22', 'k23', 'k24'};
            metList = {'*A', '*B', '*C', '*D', '*I', '*P1', '*P2', '*Q', '*R'};
            numTerm = {'+k09', '-k10.*P1'; '+k21', '-k22.*P2'};
            prodNum = [5 6; 10 11];
            reactionName = 'r_r31';
            promiscuousRxnI = 1;
            
            buildReaction(state,rateList,metList,numTerm,prodNum,reactionName, promiscuousRxnI)
          
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildReactionPromiscuous1.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);   
        end
        
                
        function testBuildReactionPromiscuous2(testCase)
             state = {'E1= k02*k05*k08*k11+k02*k05*k08*k10+k02*k05*k11*k09+k02*k08*k11*k04+k02*k08*k04*k10+k02*k11*k04*k09+k05*k08*k11*k03+k05*k08*k03*k10+k05*k11*k03*k09;', ...
                      'E2= k01.*A*k05*k08*k11+k01.*A*k05*k08*k10+k01.*A*k05*k11*k09+k01.*A*k04*k08*k11+k01.*A*k04*k08*k10+k01.*A*k04*k11*k09+k04*k06.*P1*k08*k11+k04*k06.*P1*k08*k10+k04*k06.*P1*k11*k09;', ...
                      'E3= k06.*P1*k02*k08*k11+k06.*P1*k02*k08*k10+k06.*P1*k02*k11*k09+k03*k01.*A*k08*k11+k03*k01.*A*k08*k10+k03*k01.*A*k11*k09+k06.*P1*k03*k08*k11+k06.*P1*k03*k08*k10+k06.*P1*k03*k11*k09;', ...
                      'E4= k07.*B*k02*k05*k11+k07.*B*k10*k02*k05+k10*k12.*P2*k02*k05+k07.*B*k02*k11*k04+k07.*B*k10*k02*k04+k10*k12.*P2*k02*k04+k07.*B*k05*k11*k03+k07.*B*k10*k05*k03+k10*k12.*P2*k05*k03;', ...
                      'E5= k12.*P2*k02*k05*k08+k09*k07.*B*k02*k05+k12.*P2*k09*k02*k05+k12.*P2*k02*k08*k04+k09*k07.*B*k02*k04+k12.*P2*k09*k02*k04+k12.*P2*k05*k08*k03+k09*k07.*B*k05*k03+k12.*P2*k09*k05*k03;'};

            rateList = {'k01', 'k02', 'k03', 'k04', 'k05', 'k06', 'k07', 'k08', ...
                        'k09', 'k10', 'k11', 'k12'};
            metList = {'*A', '*B', '*P1', '*P2'};
            numTerm = {'+k05', '-k06.*P1'; '+k11', '-k12.*P2'};
            prodNum = [3 1; 5 1];
            reactionName = 'r_r51';
            promiscuousRxnI = 2;
            
            buildReaction(state,rateList,metList,numTerm,prodNum,reactionName, promiscuousRxnI)
          
            filepath = fullfile(testCase.tempReactionsFolder, [reactionName, '.m']);
            res = textread(filepath,'%s');
            
            filepath = fullfile(testCase.currentPath{1}, 'testFiles', 'trueResBuildReactionPromiscuous2.txt');
            trueRes = textread(filepath,'%s');
            
            testCase.verifyEqual(res, trueRes);              
        end

    end
end
