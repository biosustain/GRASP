function [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints)
% Applies Thermodynamic Flux Balance Analysis (TMFA) to find
% thermodynamically feasible ranges for reactions fluxes, metabolite
% concentrations, and reactions Gibbs energies.
%
% It is based on https://arxiv.org/pdf/1803.04999.pdf
%
% The standard Gibbs energies are converted into metabolite
% formation energies, so that the dG values of reactions 
% involved in a loop sum to zero.
% Note that those formation energies are compatible with the
% given standard Gibbs energies, but are not necessarily
% similar to formation energies obtained from e.g. eQuilibrator.
%
% Note: The value for K is now decreased if the initial point to use for
%       sampling gibbs energies and fluxes is thermodynamically infeasible
%       (Gibbs energies and reaction fluxes are not compatible).
%       This happens sometimes because of numerical precision issues.
%       For instance, intlinprog considers anything lower than 
%       1e-6 to be the same as zero. Thus, it is possible to have 
%       a constraint like -dGr -K*e <= -1e-6, if dGr = -10, K = 1e10,
%       and e= 1e-7 the constraint is satisfied, even though it wouldn't
%       be if e=0. In this case e should be 1, which is not the case 
%       due to numerical imprecision. Thus, in this case, K would be set 
%       to 1e7, so that K*e >= 10^0 is true even when e=1e-7, forcing e 
%       to be 1 when it should be 1, e.g. to satisfy -dGr -K*e <= -1e-6,
%       with dGr = -10, K = 1e7. 
%
%
% USAGE:
%
%    [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints)
%
% INPUT:
%    ensemble (struct):             model ensemble, see buildEnsemble for fields description
%    DGf_std_min (double vector):   minimum values for metabolite formation energies
%    DGf_std_max (double vector):   maximum values for metabolite formation energies
%    vmin (double vector):          minimum values for reaction fluxes
%    vmax (double vector):          maximum values for reaction fluxes
%    xmin (double vector):	        minimum values for metabolite concentrations
%    xmax (double vector):          maximu values for metabolite concentrations
%    ineqConstraints (double):	    inequality constraints
%
% OUTPUT:
%    v_range (double matrix):        range of feasible reaction fluxes
%    DGr_range (double matrix):	     range of feasible Gibbs energies
%    DGf_std_range (double matrix):	 range of feasible formation energies
%    lnx_range (double matrix):      range of feasible metabolite concentrations (ln)
%    x0 (double matrix):             initial point to be used for the sampling with hit-and-run method
%
% .. Authors:
%       - Pedro Saa                 2016 original code
%       - Marta Matos	            2018 added error messanges and findProblematicReactions
%       - Pedro Saa and Marta Matos	2020  modified the TMFA MILP to use formation energies


Sthermo   = ensemble.Sthermo;
Sflux     = ensemble.Sflux;

% First we estimate feasible ranges for DGf_std from DGr_std_max. This is
% done solving various LP problems to ensure consistency
[m,n]       = size(Sthermo);
[~,nflux]   = size(Sflux);
M           = 1e3;
DGf_std_min = zeros(m,1);
DGf_std_max = zeros(m,1);

if strcmp(ensemble.LPSolver, 'gurobi')
    
    gurobiModel.A   = sparse([eye(n),-Sthermo']);
    gurobiModel.rhs = zeros(n,1);
    gurobiModel.lb  = [DGr_std_min(:);-M*ones(m,1)];    % Only negative DGf_std
    gurobiModel.ub  = [DGr_std_max(:);zeros(m,1)];
    gurobiModel.obj = zeros(m+n,1);
    gurobiModel.sense = '=';            % Constraint type definition
    gurobiModel.vtype = 'C';            % Variable type definition
    
    % Define optimization parameters
    params.outputflag     = 0;
    params.FeasibilityTol = 1e-9;
    params.OptimalityTol  = 1e-9;
       
    % Find feasible ranges for each DGf_std
    for ix = 1:m
        gurobiModel.obj(n+ix)  = 1;
        gurobiModel.modelsense = 'min';
        solmin = gurobi(gurobiModel,params);
        DGf_std_min(ix) = solmin.objval;
        
        gurobiModel.modelsense = 'max';
        solmax = gurobi(gurobiModel,params);
        DGf_std_max(ix) = solmax.objval;
        gurobiModel.obj(n+ix)  = 0;
    end
    
elseif strcmp(ensemble.LPSolver, 'linprog')
    
    A    = [];
    b    = [];
    Aeq  = [eye(n),-Sthermo'];
    beq  = zeros(n,1);
    lb   = [DGr_std_min(:);-M*ones(m,1)];
    ub   = [DGr_std_max(:);zeros(m,1)];
    x0   = [];
    fobj = zeros(m+n,1);
    
    options =  optimoptions(@linprog, ...
        'OptimalityTolerance', 1e-9, ...
        'ConstraintTolerance', 1e-9, ...
        'Display', 'off');
    
    % Find feasible ranges for each DGf_std
    for ix = 1:m
        fobj(n+ix) = 1;
        [~,fval]   = linprog(fobj, A, b, Aeq, beq, lb, ub, x0, options);
        DGf_std_min(ix) = fval;
        fobj(n+ix) = -1;
        [~,fval]   = linprog(fobj, A, b, Aeq, beq, lb, ub, x0, options);
        DGf_std_max(ix) = -fval;
        fobj(n+ix) = 1;
    end
    
end


% Build the adapted TMFA problem
K       = 1e12;
delta   = 1e-6;
RT      = 8.314*298.15/1e3;  % [kJ/mol]
M       = 1e5;
tol     = 1e-10;

fluxdGInconsistent = 1;

while fluxdGInconsistent
    
    % Define bounds
    [m,n]     = size(Sthermo);
    [~,nflux] = size(Sflux);
    model.lb  = [vmin;-M*ones(n,1);DGf_std_min;log(xmin);zeros(n,1)];
    model.ub  = [vmax;M*ones(n,1);DGf_std_max;log(xmax);ones(n,1)];
    Vblock    = eye(nflux);
    Vblock    = Vblock(ensemble.idxNotExch,:);


    model.Aeq = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m)]);                        % Sflux*v = 0

    % Define problem matrix
    model.A  = sparse([zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n);...     % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
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
    model.f = zeros(size(model.A,2),1);

    % Constraints sense and rhs
    model.beq = zeros(size(Sflux,1),1);
    model.b = [tol*ones(n,1);tol*ones(n,1);K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];

    if ~isempty(ineqConstraints)
        model.b = [model.b;ineqConstraints(:,end)];
    end

    if strcmp(ensemble.LPSolver, 'gurobi')

        gurobiModel.A = [model.Aeq; model.A];
        gurobiModel.rhs = [model.beq; model.b];
        gurobiModel.lb = model.lb;
        gurobiModel.ub = model.ub;
        gurobiModel.obj = model.f;

        gurobiModel.sense = blanks(numel(gurobiModel.rhs));
        gurobiModel.sense(1:size(Sflux,1)) = '=';
        gurobiModel.sense(size(Sflux,1)+1:end) = '<';

        % Variable type definition
        gurobiModel.vtype = blanks(numel(gurobiModel.obj));
        gurobiModel.vtype(1:(nflux+n+2*m)) = 'C';
        gurobiModel.vtype((nflux+n+2*m+1):end) = 'B';

        % Define optimization parameters
        params.outputflag     = 0;
        params.FeasibilityTol = 1e-6;
        params.IntFeasTol     = 1e-6;
        params.MIPGap         = 1e-6;
        params.MIPGapAbs      = 1e-12;
        params.OptimalityTol  = 1e-9;


        % Check the feasibility of the problem
        sol = gurobi(gurobiModel,params);

    elseif strcmp(ensemble.LPSolver, 'linprog')
        options =  optimoptions(@intlinprog, ...
                                'AbsoluteGapTolerance', 1e-12, ...               % equivalent to MIPGapAbs
                                'RelativeGapTolerance', 1e-6, ....               % equivalent to MIPGap
                                'LPOptimalityTolerance', 1e-9, ...               % equivalent to OptimalityTol
                                'ConstraintTolerance', 1e-6, ...                 % equivalent to FeasibilityTol
                                'IntegerTolerance', 1e-6, ...                    % equivalent to IntFeasTol
                                'Display', 'off');

        model.intcon = [(nflux+n+2*m+1):numel(model.f)];

        [x,fval] = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
    end


    if (strcmp(ensemble.LPSolver, 'gurobi') && strcmp(sol.status,'INFEASIBLE')) || ...
            (strcmp(ensemble.LPSolver, 'linprog') && isempty(fval))

        disp('TMFA problem is infeasible.')
        if strcmp(ensemble.LPSolver, 'gurobi')
            modelWithIssues = gurobiModel;
        else
            modelWithIssues = model;
        end
        
        findIssuesWithTMFA(ensemble,modelWithIssues,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol);
    end


    % Run improved TMFA
    v_range   = zeros(nflux,2);
    DGr_range = zeros(n,2);
    lnx_range = zeros(m,2);
    DGf_std_range = zeros(m,2);

    xprev_min = zeros(nflux+n+2*m,1);
    xprev_max = zeros(nflux+n+2*m,1);

    for ix = 1:nflux+n+2*m

         if strcmp(ensemble.LPSolver, 'gurobi')

            gurobiModel.obj(ix)    = 1;
            gurobiModel.modelsense = 'min';
            solmin                 = gurobi(gurobiModel,params);
            if ~strcmp(solmin.status,'OPTIMAL')
                  error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
            end

            gurobiModel.modelsense = 'max';
            solmax                 = gurobi(gurobiModel,params);
            if ~strcmp(solmax.status,'OPTIMAL')
                  error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
            end
            
            gurobiModel.obj(ix) = 0;
            solmaxX             = solmax.x;
            solminX             = solmin.x;
            solmax              = solmax.objval;
            solmin              = solmin.objval;

        elseif strcmp(ensemble.LPSolver, 'linprog')

            model.f(ix)      = 1;
            [solminX, solmin]   = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
            if  isempty(solmin)
                  error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
            end

            model.f(ix)     = -1;
            [solmaxX, solmax]  = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
            solmax          = -solmax;
            if  isempty(solmax)
                  error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
            end

            model.f(ix)     = 0;

        end  

        xprev_min(:,ix) = solminX(1:nflux+n+2*m);                            % save previous min solution for later
        xprev_max(:,ix) = solmaxX(1:nflux+n+2*m);                            % save previous max solution for later

        if (ix<=nflux)
            v_range(ix,:) = [solmin, solmax];                    % [umol or mmol/gdcw/h]
        elseif (ix<=nflux+n)
            DGr_range(ix-nflux,:) = [solmin, solmax];            % [kJ/mol]
        elseif (ix<=nflux+n+m)
            DGf_std_range(ix-nflux-n,:) = [solmin, solmax];      % [kJ/mol]
        else
            lnx_range(ix-nflux-n-m,:)   = [solmin, solmax];      % log([M])
        end

    end


    % Return initial point
    x0 = mean([xprev_min,xprev_max],2);

    Nint = null(Sthermo,'r');
    loopCondition = Nint' * x0(nflux+1:nflux+n);

    % Check that the initial point is feasible
    model.A = [model.Aeq; model.A];
    Acheck = full(model.A(1:(size(Sflux,1)+2*n),1:(nflux+n+2*m)));              % Generate checking matrix

    if all(abs(Acheck*x0)<1e-6)
         
        %Check that the initial point fluxes and dG values have opposite signs
         if sum(sign(x0(ensemble.idxNotExch)) + sign(x0(ensemble.idxNotExch+nflux))) == 0 

            disp('The fluxes and Gibbs energies are consistent.');
            fluxdGInconsistent = 0;
            
            if K < 1e12
                disp(['The K value in the TMFA problem was modified and set to 10^', num2str(log10(K))]);
            end

         else
            n = log10(K);
            K = 10^(n-1);
            if K >= max(abs(v_range), [], 'all') && K >= max(abs(DGr_range), [], 'all')
                continue
            else
               error('The fluxes and Gibbs energies are not consistent.');
            end
         end
         
        if sum(sign(xprev_min(1:nflux)) - sign(xprev_max(1:nflux))) == 0    
            disp('Flux directions are consistent.');
        else 
            error('The flux directions are not consistent. Please make sure that both the lower and upper bound of the flux ranges (fluxMean - 2*fluxStd and fluxMean + 2*fluxStd, respectively) are either positive or negative.');
        end
        
        nonZeroIndMin = find(xprev_min(nflux+1:nflux+n) ~= 0);
        nonZeroIndMax = find(xprev_max(nflux+1:nflux+n) ~= 0);
        nonZeroIntersec = intersect(nonZeroIndMin, nonZeroIndMax);

        if sum(sign(xprev_min(nonZeroIntersec)) - sign(xprev_max(nonZeroIntersec))) == 0    
            disp('Gibbs energy directions are consistent.');
        else 
            error('The mininum and maximum values of Gibbs energies used to calculate the initial point for dG and flux sampling are not consistent.');
        end


        if ~isempty(Nint)
            if max(abs(loopCondition)) <= 1e-6
                disp('The initial point obtained from TMFA is feasible and valid for starting the sampler.');
            else
                error('The initial point obtained from TMFA is not thermodynamically feasible.')
            end
        else
            disp('The initial point obtained from TMFA is feasible and valid for starting the sampler.');
        end            

    else
        error('The initial point obtained from TMFA is not thermodynamically feasible.')
    end
end

end
