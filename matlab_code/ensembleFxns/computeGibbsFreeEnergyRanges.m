function [DGr_rng,xrng,vrng] = computeGibbsFreeEnergyRanges(Sflux,Sthermo,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints, rxnNames, solver)
% Applies Thermodynamic Flux Balance Analysis (TMFA) to find
% thermodynamically feasible ranges for reactions fluxes, metabolite
% concentrations, and reactions Gibbs energies.
%
% It is based on https://arxiv.org/pdf/1803.04999.pdf
%
%
% USAGE:
%
%    [DGr_rng, xrng, vrng] = computeGibbsFreeEnergyRanges(Sflux, Sthermo, DGr_std_min, DGr_std_max, vmin, vmax, xmin, xmax, idxNotExch, ineqConstraints, rxnNames)
%
% INPUT:
%    Sflux (int matrix):          stoichiometric matrix used for flux calculations
%    Sthermo (int matrix):        stoichiometric matrix used for  thermodynamic calculations
%    DGr_std_min (double vector): lower bound for standard Gibbs energies
%    DGr_std_max (double):		  upper bound for standard Gibbs energies
%    vmin (double):               lower bound for reactions fluxes
%    vmax (double):               upper bound for reactions fluxes
%    xmin (double):               lower bound for metabolite concentrations
%    xmax (double):               upper bound for metabolite concentrations
%    idxNotExch (int vector):     IDs for reactions that are not exchange reactions
%    ineqConstraints (double):	  inequality constraints
%    rxnNames (char vector):      reaction names
%    solver (char vector):        specifies which solver to use to solve the MILP: 'gurobi' or 'linprog'
%
% OUTPUT:
%    DGr_rng (double matrix):	range of feasible Gibbs energies
%    xrng (double matrix):      range of feasible metabolite concentrations
%    vrng (double matrix):      range of feasible reaction fluxes
%
% .. Authors:
%       - Pedro Saa     2016 original code
%       - Marta Matos	2018 added error messanges and 
%                       findProblematicReactions

% Build the adapted TMFA problem
K       = 1e5;
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
model.Aeq = sparse([Sflux,zeros(size(Sflux,1),2*n+m)]);                   % Sflux*v = 0

model.A  = sparse([zeros(n,nflux),eye(n),-RT*Sthermo',zeros(n);...        % DGr - RT*Sthermo*ln(x) <= DGr_std_max
                   zeros(n,nflux),-eye(n),RT*Sthermo',zeros(n);...        % -DGr + RT*Sthermo*ln(x) <= -DGr_std_min
                   -Vblock,zeros(n,n+m),K*eye(n);...                      % -v + K*e <= K
                   Vblock,zeros(n,n+m),-K*eye(n);...                      % v - K*e <= 0
                   zeros(n,nflux),eye(n),zeros(n,m),K*eye(n);...          % DGr + K*e <= K - delta
                   zeros(n,nflux),-eye(n),zeros(n,m),-K*eye(n)]);         % -DGr - K*e <= -delta
         
if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.f = zeros(size(model.A,2),1);

% Constraints sense and rhs
model.beq = zeros(size(Sflux,1),1);
model.b = [DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];

if ~isempty(ineqConstraints)
    model.b = [model.b;ineqConstraints(:,end)];
end

if strcmp(solver, 'gurobi')
    
    gurobiModel.A = [model.Aeq; model.A];
    gurobiModel.rhs = [model.beq; model.b];
    gurobiModel.lb = model.lb;
    gurobiModel.ub = model.ub;
    gurobiModel.obj = model.f;
    
    gurobiModel.sense = blanks(numel(gurobiModel.rhs));
    for ix = 1:numel(gurobiModel.rhs)
        if (ix<=size(Sflux,1))
            gurobiModel.sense(ix) = '=';
        else
            gurobiModel.sense(ix) = '<';
        end
    end

    % Variable type definition
    gurobiModel.vtype = blanks(numel(gurobiModel.obj));
    for ix = 1:numel(gurobiModel.obj)
        if (ix<=nflux+n+m)
            gurobiModel.vtype(ix) = 'C';
        else
            gurobiModel.vtype(ix) = 'B';
        end
    end
    
    % Define optimization parameters
    params.outputflag = 0;
    params.OptimalityTol = 1e-6;
    params.FeasibilityTol = 1e-6;
    params.IntFeasTol = 1e-5;


    % Check the feasibility of the problem
    sol = gurobi(gurobiModel,params);

elseif strcmp(solver, 'linprog')
    options =  optimoptions(@intlinprog, 'LPOptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
                            'IntegerTolerance', 1e-6, 'Display', 'off');
                        
    model.intcon = [(nflux+n+m+1):numel(model.f)];
    
    [x,fval] = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
end


if (strcmp(solver, 'gurobi')&& strcmp(sol.status,'INFEASIBLE')) || ...
        (strcmp(solver, 'linprog') && isempty(fval))
    
    [row_list, dg_list] = findProblematicReactions(gurobiModel, params, model, options, DGr_std_min, DGr_std_max, K, delta, n, Sflux, ineqConstraints, rxnNames, solver);
   
    if isempty(row_list)
        error('The TMFA problem is infeasible but it is not due to incompatible fluxes and Gibbs energies. It is likely to be an issue with fluxes. Make sure v - K <= 0, where K = 1e5.');
    else
        error(strcat('The TMFA problem is infeasible. Verify that the standard Gibbs free energy and metabolite concentration values are valid/correct. Reactions', strjoin(row_list, ', '), ' with standard Gibbs energies ', mat2str(dg_list), ' seem to be the problem. Note that the problem can also be with the reaction fluxes. Bottom line, for each reaction, fluxes and Gibbs energies need to agree.'));
    end    
end

% Run improved TMFA
vrng    = zeros(nflux,2);
DGr_rng = zeros(n,2);
xrng    = zeros(m,2);
for ix = 1:nflux+n+m
    
    if strcmp(solver, 'gurobi')

        gurobiModel.obj(ix)    = 1;
        gurobiModel.modelsense = 'min';
        solmin           = gurobi(gurobiModel,params);
        if ~strcmp(solmin.status,'OPTIMAL')
              error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
        end
        
        gurobiModel.modelsense = 'max';
        solmax           = gurobi(gurobiModel,params);
        
        gurobiModel.obj(ix) = 0;
        solmax           = solmax.objval;
        solmin           = solmin.objval;
    
    elseif strcmp(solver, 'linprog')
        
        model.f(ix)            = 1;
        [x, solmin]      = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
        if  isempty(solmin)
              error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
        end
        
        model.f(ix)           = -1;
        [x, solmax]     = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
        solmax          = -solmax;
        model.f(ix)           = 0;
        
    end
    if (ix<=nflux)
        vrng(ix,:)          = [solmin,solmax];      % [flux units as give in measRates sheet]
    elseif (ix<=nflux+n)
        DGr_rng(ix-nflux,:) = [solmin,solmax];      % [kJ/mol]
    else
        xrng(ix-nflux-n,:)  = [solmin,solmax];      % exp([M])
    end
    
end
xmax = max(xrng(:));  % robust calculation of exp(log(x))
xrng = xrng - xmax;
xrng = exp(xrng);
xrng = xrng*exp(xmax);
