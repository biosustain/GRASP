classdef computeGibbsFreeEnergyRangesTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
        end
    end
    
 
    methods (Test)
        function testComputeGibbsFreeEnergyRanges1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble =  load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1');
            xlsxFile            = [xlsxFile,'.xlsx'];                               % add extension to the file
            xDG_std             = xlsread(xlsxFile,'thermoRxns');                   % load thermodynamic data (rxns)
            xMetsThermo         = xlsread(xlsxFile,'thermoMets');                   % load thermodynamic data (mets)
            ineqConstraints     = xlsread(xlsxFile,'thermo_ineq_constraints');      % load ineq. thermodynamic constraints

            Sflux       = ensemble.S(ensemble.metsBalanced,:);
            idxNotExch  = find(~ismember(1:numel(ensemble.rxns),ensemble.exchRxns));
            ensemble.Sthermo     = ensemble.S(:,idxNotExch);
            DGr_std     = xDG_std(idxNotExch,:);                                                        % Use only reactions with known thermodynamics
            vmin        = ensemble.fluxRef - 2*ensemble.fluxRefStd;
            vmax        = ensemble.fluxRef + 2*ensemble.fluxRefStd;
            xmin        = xMetsThermo(:,1);
            xmax        = xMetsThermo(:,2);
            DGr_std_min = DGr_std(:,1);
            DGr_std_max = DGr_std(:,2);
            
            [DGr_rng,xrng,vrng] = computeGibbsFreeEnergyRanges(Sflux,ensemble.Sthermo,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints', ensemble.rxnNames);
           
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResComputeGibbsFreeEnergyRanges1'));
            trueResDGr = trueRes.DGr_rng;
            trueResMets = trueRes.xrng;
            trueResFluxes = trueRes.vrng;
                   
            testCase.verifyEqual(trueResDGr, DGr_rng);
            testCase.verifyEqual(trueResMets, xrng);
            testCase.verifyEqual(trueResFluxes, vrng);
        end
       
    end
end

