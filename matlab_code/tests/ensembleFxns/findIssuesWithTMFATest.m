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

            DGf_std_max = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_max;
            DGf_std_min = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_min;
            for ix = 1:numel(DGf_std_max)
                if (DGf_std_min(ix)>DGf_std_max(ix))
                    DGf_std_temp    = DGf_std_max(ix);
                    DGf_std_max(ix) = DGf_std_min(ix);
                    DGf_std_min(ix) = DGf_std_temp;
                end
            end


            % Build the adapted TMFA problem
            K       = 1e12;
            delta   = 1e-6;
            RT      = 8.314*298.15/1e3;  % [kJ/mol]
            M       = 1e5;
            tol     = 1e-10;

            % Define bounds
            [m,n]     = size(Sthermo);
            [~,nflux] = size(Sflux);
            model.lb  = [vmin;-M*ones(n,1);DGf_std_min;log(xmin);zeros(n,1)];
            model.ub  = [vmax;M*ones(n,1);DGf_std_max;log(xmax);ones(n,1)];
            Vblock    = eye(nflux);
            Vblock    = Vblock(ensemble.idxNotExch,:);

            % Define problem matrix
            model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m);...         % Sflux*v = 0
                zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n);...     % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
                zeros(n,nflux),-eye(n),Sthermo',RT*Sthermo',zeros(n);...      % -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) <= tol
                -Vblock,zeros(n,n+2*m),K*eye(n);...                           % -v + K*e <= K
                Vblock,zeros(n,n+2*m),-K*eye(n);...                           % v - K*e <= 0
                zeros(n,nflux),eye(n),zeros(n,2*m),K*eye(n);...               % DGr + K*e <= K - delta
                zeros(n,nflux),-eye(n),zeros(n,2*m),-K*eye(n)]);              % -DGr - K*e <= -delta

            if ~isempty(ineqConstraints)
                p = size(ineqConstraints,1);
                model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
            end

            % Objective function
            model.obj = zeros(size(model.A,2),1);

            % Constraints sense and rhs
            model.rhs = [zeros(size(Sflux,1),1);tol*ones(n,1);tol*ones(n,1);K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];
            if ~isempty(ineqConstraints)
                model.rhs = [model.rhs;ineqConstraints(:,end)];
            end

            model.sense = blanks(numel(model.rhs));
            model.sense(1:size(Sflux,1)) = '=';
            model.sense(size(Sflux,1)+1:end) = '<';

            % Variable type definition
            model.vtype = blanks(numel(model.obj));
            model.vtype(1:nflux+n+2*m) = 'C';
            model.vtype(nflux+n+2*m+1:end) = 'B';

            testCase.verifyError(@()findIssuesWithTMFA(ensemble,model,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol), '');   
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

            DGf_std_max = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_max;
            DGf_std_min = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_min;
            for ix = 1:numel(DGf_std_max)
                if (DGf_std_min(ix)>DGf_std_max(ix))
                    DGf_std_temp    = DGf_std_max(ix);
                    DGf_std_max(ix) = DGf_std_min(ix);
                    DGf_std_min(ix) = DGf_std_temp;
                end
            end


            % Build the adapted TMFA problem
            K       = 1e12;
            delta   = 1e-6;
            RT      = 8.314*298.15/1e3;  % [kJ/mol]
            M       = 1e5;
            tol     = 1e-10;
            
            ensemble.LPSolver = 'linprog';
            model = [];
            
            testCase.verifyError(@()findIssuesWithTMFA(ensemble,model,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol), '');   
        end
       
    end
end