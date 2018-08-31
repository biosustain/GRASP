function [DrG_ranges,logX_ranges] = calculateGibbsEnergyRanges(Sfull,DrG_std,metRanges,dirVector,ineqThermoConst)
%--------------------------------------------------------------------------
% Performs thermodynamic variability analysis
% Inputs:  Stoichiometric matrix (without exchs), DrG_std (ranges), metRanges, directionality vector
%
% Outputs: Gibbs free energy ranges (kJ/mol)
%
%------------------------Pedro Saa 2016------------------------------------
% Optimization parameters
K  = 1e2;                        % large constant
RT = 8.314*298.15/1e3;           % gas constant*temp (assumed 25 C temp.)
params.outputflag     = 0;       % gurobi params
params.OptimalityTol  = 1e-9;
params.FeasibilityTol = 1e-9;

% Formulate TVA problem (LP formulation)
[m,n]           = size(Sfull);
model.obj       = zeros(2*n+m,1);
model.A         = sparse([eye(n),-eye(n),-RT*Sfull']);
if isempty(ineqThermoConst)
    model.rhs   = zeros(n,1);
    model.sense = '=';
else
    p           = size(ineqThermoConst,1);    
    model.rhs   = zeros(n,1);
    model.rhs   = [model.rhs;ineqThermoConst(:,end)];
    model.A     = sparse([model.A;zeros(p,n),zeros(p,n),ineqThermoConst(:,1:end-1)]);
    model.sense = blanks(n+p);
    for ix = 1:n+p
        if (ix <= n)
            model.sense(ix) = '=';
        else
            model.sense(ix) = '<';
        end
    end
end

model.vtype      = 'C';
model.modelsense = 'min';
model.lb         = [-dirVector/K-K*(1+dirVector);DrG_std(:,1);log(metRanges(:,1))];
model.ub         = [-dirVector/K+K*(1-dirVector);DrG_std(:,2);log(metRanges(:,2))];
DrG_ranges       = zeros(n,2);
f                = model.obj;

for ix = 1:n
    model.obj(ix) = 1;
    
    % Minimization problem
    model.modelsense = 'min';
    solmin           = gurobi(model,params);
    assert(~strcmp(solmin.status, 'INFEASIBLE'), 'TMFA problem is infeasbile, most likely there is a conflict between the Gibbs energies and the reaction fluxes');   
    DrG_ranges(ix,1) = solmin.objval;
    
    % Maximization problem
    model.modelsense = 'max';
    solmax           = gurobi(model,params);
    DrG_ranges(ix,2) = solmax.objval;
    
    % Restart obj fxn
    model.obj = f;
end

logX_ranges = zeros(m,2);
for ix = 2*n+1:2*n+m
    model.obj(ix) = 1;
    
    % Minimization problem
    model.modelsense = 'min';
    solmin           = gurobi(model,params);
    logX_ranges(ix,1) = solmin.objval;
    
    % Maximization problem
    model.modelsense = 'max';
    solmax           = gurobi(model,params);
    logX_ranges(ix,2) = solmax.objval;
    
    % Restart obj fxn
    model.obj = f;
end

logX_ranges((logX_ranges(:,1)==0)&(logX_ranges(:,2)==0),:) = [];


% Solve for maximum and minimum ratio of fld_red/fld_ox
%model.obj([35;36]) = [-1;1];

% Minimization problem
%model.modelsense  = 'min';
%solmin            = gurobi(model,params);
%minimum_logratio_fld = solmin.objval

% Maximization problem
%model.modelsense = 'max';
%solmax           = gurobi(model,params);
%maximum_logratio_fld = solmax.objval
