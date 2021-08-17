classdef findIssuesWithTMFATest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testFindIssuesWithTMFAGurobi(testCase)
            
            seed = 1;
            rng(seed)
           
            ensemble =  load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model2_test_gibbs.mat'));
            ensemble = ensemble.ensemble;

            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model2_test_gibbs');
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
            
            Sthermo   = ensemble.Sthermo;
            Sflux     = ensemble.Sflux;

            % Build the adapted TMFA problem
            K       = 1e12;
            delta   = 1e-6;
            RT      = 8.314*298.15/1e3;  % [kJ/mol]
            M       = 1e5;
            tol     = 1e-10;

            gurobiModel= [];
                        
            options =  optimoptions(@intlinprog, ...
                                    'AbsoluteGapTolerance', 1e-12, ...               % equivalent to MIPGapAbs
                                    'RelativeGapTolerance', 1e-6, ....               % equivalent to MIPGap
                                    'LPOptimalityTolerance', 1e-9, ...               % equivalent to OptimalityTol
                                    'ConstraintTolerance', 1e-6, ...                 % equivalent to FeasibilityTol
                                    'IntegerTolerance', 1e-6, ...                    % equivalent to IntFeasTol
                                    'Display', 'off');
            
            params.outputflag     = 0;
            params.FeasibilityTol = 1e-6;
            params.IntFeasTol     = 1e-6;
            params.MIPGap         = 1e-6;
            params.MIPGapAbs      = 1e-12;
            params.OptimalityTol  = 1e-9;

            testCase.verifyError(@()findIssuesWithTMFA(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,gurobiModel,options,params), '');   
        end
        
        function testFindIssuesWithTMFALinprog(testCase)
            
            ensemble =  load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model2_test_gibbs.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model2_test_gibbs');
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
            
            Sthermo   = ensemble.Sthermo;

            % Build the adapted TMFA problem
            K       = 1e12;
            delta   = 1e-6;
            RT      = 8.314*298.15/1e3;  % [kJ/mol]
            M       = 1e5;
            tol     = 1e-10;
            
            ensemble.LPSolver = 'linprog';
            
            gurobiModel = [];
            
            options =  optimoptions(@intlinprog, ...
                                    'AbsoluteGapTolerance', 1e-12, ...               % equivalent to MIPGapAbs
                                    'RelativeGapTolerance', 1e-6, ....               % equivalent to MIPGap
                                    'LPOptimalityTolerance', 1e-9, ...               % equivalent to OptimalityTol
                                    'ConstraintTolerance', 1e-6, ...                 % equivalent to FeasibilityTol
                                    'IntegerTolerance', 1e-6, ...                    % equivalent to IntFeasTol
                                    'Display', 'off');
            
            params.outputflag     = 0;
            params.FeasibilityTol = 1e-6;
            params.IntFeasTol     = 1e-6;
            params.MIPGap         = 1e-6;
            params.MIPGapAbs      = 1e-12;
            params.OptimalityTol  = 1e-9;
            
            testCase.verifyError(@()findIssuesWithTMFA(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,gurobiModel,options,params), '');  
        end
       
    end
end