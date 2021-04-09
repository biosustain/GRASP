classdef computeGibbsFreeEnergyRangesTest < matlab.unittest.TestCase

    properties
        currentPath
        relTol = 1e-4;
        absTol = 1e-4;
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testComputeGibbsFreeEnergyRangesGurobi(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1.xlsx');
            xDG_std             = xlsread(xlsxFile,'thermoRxns');                   % load thermodynamic data (rxns)
            xMetsThermo         = xlsread(xlsxFile,'thermoMets');                   % load thermodynamic data (mets)
            ineqConstraints     = xlsread(xlsxFile,'thermo_ineq_constraints');      % load ineq. thermodynamic constraints

            DGr_std        = xDG_std(ensemble.idxNotExch,:);                                                        % Use only reactions with known thermodynamics
            vmin           = ensemble.fluxRef - 2*ensemble.fluxRefStd;
            vmax           = ensemble.fluxRef + 2*ensemble.fluxRefStd;
            xmin           = xMetsThermo(:,1);
            xmax           = xMetsThermo(:,2);
            DGr_std_min    = DGr_std(:,1);
            DGr_std_max    = DGr_std(:,2);
            
            ensemble.LPSolver = 'gurobi';
            
            [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResComputeGibbsFreeEnergyRanges1'));
            trueResFluxes = trueRes.v_range;
            trueResDGr = trueRes.DGr_range;
            trueResDGf = trueRes.DGf_std_range;
            trueResLnMets = trueRes.lnx_range;
            trueResX0 = trueRes.x0;
            
            testCase.verifyThat(v_range, matlab.unittest.constraints.IsEqualTo(trueResFluxes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(DGr_range, matlab.unittest.constraints.IsEqualTo(trueResDGr, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(DGf_std_range, matlab.unittest.constraints.IsEqualTo(trueResDGf, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(lnx_range, matlab.unittest.constraints.IsEqualTo(trueResLnMets, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(x0, matlab.unittest.constraints.IsEqualTo(trueResX0, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-2) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
        end
        
        function testComputeGibbsFreeEnergyRangesLinprog(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model1.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1.xlsx');
            xDG_std             = xlsread(xlsxFile,'thermoRxns');                   % load thermodynamic data (rxns)
            xMetsThermo         = xlsread(xlsxFile,'thermoMets');                   % load thermodynamic data (mets)
            ineqConstraints     = xlsread(xlsxFile,'thermo_ineq_constraints');      % load ineq. thermodynamic constraints

            DGr_std        = xDG_std(ensemble.idxNotExch,:);                                                        % Use only reactions with known thermodynamics
            vmin           = ensemble.fluxRef - 2*ensemble.fluxRefStd;
            vmax           = ensemble.fluxRef + 2*ensemble.fluxRefStd;
            xmin           = xMetsThermo(:,1);
            xmax           = xMetsThermo(:,2);
            DGr_std_min    = DGr_std(:,1);
            DGr_std_max    = DGr_std(:,2);
            
            ensemble.LPSolver = 'linprog';
            
            [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResComputeGibbsFreeEnergyRangesLinprog'));
            trueResFluxes = trueRes.v_range;
            trueResDGr = trueRes.DGr_range;
            trueResDGf = trueRes.DGf_std_range;
            trueResLnMets = trueRes.lnx_range;
            trueResX0 = trueRes.x0;
            
            testCase.verifyThat(v_range, matlab.unittest.constraints.IsEqualTo(trueResFluxes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(DGr_range, matlab.unittest.constraints.IsEqualTo(trueResDGr, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(DGf_std_range, matlab.unittest.constraints.IsEqualTo(trueResDGf, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(lnx_range, matlab.unittest.constraints.IsEqualTo(trueResLnMets, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)))
            testCase.verifyThat(x0, matlab.unittest.constraints.IsEqualTo(trueResX0, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(0.2) | matlab.unittest.constraints.AbsoluteTolerance(0.2)))
        end
       
    end
end

