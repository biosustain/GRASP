classdef computRobustFluxesTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
    
   
    methods (Test)
        function testComputeRobustFluxes1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble =  load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            xlsxFile            = [xlsxFile,'.xlsx'];                               
            [measRates,idxMeas] = xlsread(xlsxFile,'measRates'); 

            Sflux                 = ensemble.S(ensemble.metsBalanced,:);
            idxMeas = fixVariableNames(idxMeas, 'r');
            idxMeas               = idxMeas(2:end,1);                  
            xMean = zeros(size(Sflux,2),1);
            xStd  = xMean;
            ix = 1;
            
            for jx = 1:size(idxMeas)                               
                index        = strcmp(ensemble.rxns,idxMeas(jx));
                xMean(index) = ensemble.measRates(jx,ix);
                xStd(index)  = ensemble.measRatesStd(jx,ix);
            end
            
            [vMean,vStd] = computeRobustFluxes(Sflux,xMean,xStd);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResComputeRobustFluxes1'));
            trueResVmean = trueRes.vMean;
            trueResVstd = trueRes.vStd;
            
            testCase.verifyThat(trueResVmean, matlab.unittest.constraints.IsEqualTo(vMean, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
            testCase.verifyThat(trueResVstd, matlab.unittest.constraints.IsEqualTo(vStd, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
       
    end
end
