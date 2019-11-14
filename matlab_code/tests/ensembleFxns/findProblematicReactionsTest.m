classdef findProblematicReactionsTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
        function testFindProblematicReactionsTest1(testCase)
            
            seed = 1;
            rng(seed)
            
            ensemble =  load(fullfile(testCase.currentPath{1}, 'testFiles', 'ensemble_toy_model2_test_gibbs.mat'));
            ensemble = ensemble.ensemble;
            
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model2_test_gibbs');
            xlsxFile            = [xlsxFile,'.xlsx'];                               % add extension to the file
            xDG_std             = xlsread(xlsxFile,'thermoRxns');                   % load thermodynamic data (rxns)
            xMetsThermo         = xlsread(xlsxFile,'thermoMets');                   % load thermodynamic data (mets)
            ineqConstraints     = xlsread(xlsxFile,'thermo_ineq_constraints');      % load ineq. thermodynamic constraints
            
            ineqConstraints = ineqConstraints';
            Sflux       = ensemble.S(ensemble.metsBalanced,:);
            idxNotExch  = find(~ismember(1:numel(ensemble.rxns),ensemble.exchRxns));
            Sthermo     = ensemble.S(:,idxNotExch);
            DGr_std     = xDG_std(idxNotExch,:);                                                        % Use only reactions with known thermodynamics
            vmin        = ensemble.fluxRef - 2*ensemble.fluxRefStd;
            vmax        = ensemble.fluxRef + 2*ensemble.fluxRefStd;
            xmin        = xMetsThermo(:,1);
            xmax        = xMetsThermo(:,2);
            DGr_std_min = DGr_std(:,1);
            DGr_std_max = DGr_std(:,2);
            
            % Build the adapted TMFA problem
            K       = 1e4;
            delta   = 1e-6;
            RT      = 8.314*298.15/1e3;  % [kJ/mol]

            % Define bounds
            [m,n]     = size(Sthermo);
            [~,nflux] = size(Sflux);
            model.lb  = [vmin;-K*ones(n,1);log(xmin);zeros(n,1)];
            model.ub  = [vmax;K*ones(n,1);log(xmax);ones(n,1)];
            Vblock    = eye(nflux);
            Vblock    = Vblock(idxNotExch,:);

            % Define problem matrix
            model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+m);...           % Sflux*v = 0
                zeros(n,nflux),eye(n),-RT*Sthermo',zeros(n);...               % DGr - RT*Sthermo*ln(x) <= DGr_std_max
                zeros(n,nflux),-eye(n),RT*Sthermo',zeros(n);...               % -DGr + RT*Sthermo*ln(x) <= -DGr_std_min
                -Vblock,zeros(n,n+m),K*eye(n);...                             % -v + K*e <= K
                Vblock,zeros(n,n+m),-K*eye(n);...                             % v - K*e <= 0
                zeros(n,nflux),eye(n),zeros(n,m),K*eye(n);...                 % DGr + K*e <= K - delta
                zeros(n,nflux),-eye(n),zeros(n,m),-K*eye(n)]);                % -DGr - K*e <= -delta
            if ~isempty(ineqConstraints)
                p = size(ineqConstraints,1);
                model.A = sparse([model.A;zeros(p,nflux+n),ineqConstraints(:,1:end-1),zeros(p,n)]);
            end

            % Objective function
            model.obj = zeros(size(model.A,2),1);

            % Constraints sense and rhs
            model.rhs = [zeros(size(Sflux,1),1);DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];
            if ~isempty(ineqConstraints)
                model.rhs = [model.rhs;ineqConstraints(:,end)];
            end
            model.sense = blanks(numel(model.rhs));
            for ix = 1:numel(model.rhs)
                if (ix<=size(Sflux,1))
                    model.sense(ix) = '=';
                else
                    model.sense(ix) = '<';
                end
            end

            % Variable type definition
            model.vtype = blanks(numel(model.obj));
            for ix = 1:numel(model.obj)
                if (ix<=nflux+n+m)
                    model.vtype(ix) = 'C';
                else
                    model.vtype(ix) = 'B';
                end
            end

            % Define optimization parameters
            params.outputflag = 0;

            % Check the feasibility of the problem
            sol = gurobi(model,params);
            
            [rowList, dGList] = findProblematicReactions(model,params,DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,ensemble.rxns);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResFindProblematicReactionsTest1'));
            trueResRowList = trueRes.rowList;
            trueResdGList = trueRes.dGList;
                   
            testCase.verifyEqual(trueResRowList, rowList);
            testCase.verifyEqual(trueResdGList, dGList);
        end
       
    end
end