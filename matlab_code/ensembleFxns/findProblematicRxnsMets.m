function [metsList,metsBoundaryList,rxnsList,rxnsBoundaryList] = findProblematicRxnsMets(ensemble,Sflux,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints,K,RT,delta,M,tol)

[metsList, metsBoundaryList] = findProblematicMetabolites(ensemble,Sflux,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints,K,RT,delta,M,tol);

[rxnsList, rxnsBoundaryList] = findProblematicReactions(ensemble,Sflux,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints,K,RT,delta,M,tol);

end


function [rxnsList, boundaryList] = findProblematicReactions(ensemble,Sflux,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints,K,RT,delta,M,tol)
% Finds the reactions whose standard Gibbs energy bounds make the 
% TMFA problem infeasible.
%
% First, it sets the bounds for the standard Gibbs energies to 
% [-100, 100] kJ/mol and solves the linear problem. Then, while the 
% solution is optimal it sets each reaction dG bounds to the original 
% values. If at any point the  problem becomes infeasible, it sets 
% the reaction dG bounds to [-100, 100] kJ/mol again.
% In the end it returns a list of all reactions whose dG bounds are 
% [-100, 100] kJ/mol, i.e. the reactions that cause the linear problem
% solution to be infeasible when their dG bounds are the original. 
%
%
% USAGE:
%
%    [rowList, dGList] = findProblematicReactions(model, params, DGr_std_min, DGr_std_max, K, delta, n, Sflux, ineqConstraints, sol, rxnNames)
%
% INPUT:
%    gurobiModel (struct):          gurobi LP model
%    params (struct):               parameters for gurobi MILP model
%    model (struct):                model to be used with intlinprog
%    options (struct):              parameters for intlinprog MILP model
%    DGr_std_min (double vector):   lower bound for standard Gibbs energies
%    DGr_std_max (double vector):	upper bound for standard Gibbs energies
%    K (double):                    arbitrarily large number for MILP model
%    delta (double):                tolerance for MILP model
%    n (double):                    number of reactions in the model
%    Sflux (double):                stoichiometric matrix used for flux calculations
%    ineqConstraints (double):	    inequality constraints
%    rxnNames (char vector):        reaction names
%    solver (char vector):          specifies which solver to use to solve the MILP: 'gurobi' or 'linprog'
%
% OUTPUT:
%    rowList (int vector):	  names of reactions that make the LP infeasible
%    dGList (double vector):  dG values of reactions that make the LP infeasible
%
% .. Authors:
%       - Marta Matos	2019 original code


% Define bounds
Sthermo   = ensemble.Sthermo;
[m,n]     = size(Sthermo);
[~,nflux] = size(Sflux);
model.lb  = [vmin;-M*ones(n,1);DGf_std_min;log(xmin);zeros(3*n,1)];
model.ub  = [vmax;M*ones(n,1);DGf_std_max;log(xmax);ones(3*n,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(idxNotExch,:);

% Define problem matrix
model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m+2*n);...                       % Sflux*v = 0
    zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n),K*eye(n),zeros(n,n);...   % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) + K*ubound <= tol
    zeros(n,nflux),-eye(n),Sthermo',RT*Sthermo',zeros(n),zeros(n,n),-K*eye(n);...   % -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) - K*lbound <= tol
    -Vblock,zeros(n,n+2*m),K*eye(n),zeros(n,2*n);...                                % -v + K*e <= K
    Vblock,zeros(n,n+2*m),-K*eye(n),zeros(n,2*n);...                                % v - K*e <= 0
    zeros(n,nflux),eye(n),zeros(n,2*m),K*eye(n),zeros(n,2*n);...                    % DGr + K*e <= K - delta
    zeros(n,nflux),-eye(n),zeros(n,2*m),-K*eye(n),zeros(n,2*n)]);                   % -DGr - K*e <= -delta

if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.obj = [zeros(size(model.A,2)-2*n,1); ones(2*n,1)];

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

% Define optimization parameters
params.outputflag    = 0;
params.IntFeasTol    = 1e-9;
params.MIPGap        = 1e-6;
params.MIPGapAbs     = 1e-12;
params.OptimalityTol = 1e-9;

% Check the feasibility of the problem
sol = gurobi(model,params);


rowList = find(sol.x(end-2*n+1:end)==1);
boundaryList = cell(numel(rowList),1);

for entry=1:numel(rowList)
    row = rowList(entry);
    
    if row > n
        row = row-n;
        rowList(entry) = row;
        boundaryList{entry} = 'lower';
    else
        boundaryList{entry} = 'upper';
    end
    
end

rxnsList = ensemble.rxns(rowList);

end



function [metList, boundaryList] = findProblematicMetabolites(ensemble,Sflux,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints,K,RT,delta,M,tol)
% Finds the reactions whose standard Gibbs energy bounds make the 
% TMFA problem infeasible.
%
% First, it sets the bounds for the standard Gibbs energies to 
% [-100, 100] kJ/mol and solves the linear problem. Then, while the 
% solution is optimal it sets each reaction dG bounds to the original 
% values. If at any point the  problem becomes infeasible, it sets 
% the reaction dG bounds to [-100, 100] kJ/mol again.
% In the end it returns a list of all reactions whose dG bounds are 
% [-100, 100] kJ/mol, i.e. the reactions that cause the linear problem
% solution to be infeasible when their dG bounds are the original. 
%
%
% USAGE:
%
%    [rowList, dGList] = findProblematicReactions(model, params, DGr_std_min, DGr_std_max, K, delta, n, Sflux, ineqConstraints, sol, rxnNames)
%
% INPUT:
%    gurobiModel (struct):          gurobi LP model
%    params (struct):               parameters for gurobi MILP model
%    model (struct):                model to be used with intlinprog
%    options (struct):              parameters for intlinprog MILP model
%    DGr_std_min (double vector):   lower bound for standard Gibbs energies
%    DGr_std_max (double vector):	upper bound for standard Gibbs energies
%    K (double):                    arbitrarily large number for MILP model
%    delta (double):                tolerance for MILP model
%    n (double):                    number of reactions in the model
%    Sflux (double):                stoichiometric matrix used for flux calculations
%    ineqConstraints (double):	    inequality constraints
%    rxnNames (char vector):        reaction names
%    solver (char vector):          specifies which solver to use to solve the MILP: 'gurobi' or 'linprog'
%
% OUTPUT:
%    rowList (int vector):	  names of reactions that make the LP infeasible
%    dGList (double vector):  dG values of reactions that make the LP infeasible
%
% .. Authors:
%       - Marta Matos	2019 original code



% Define bounds
Sthermo   = ensemble.Sthermo;
[m,n]     = size(Sthermo);
[~,nflux] = size(Sflux);
model.lb  = [vmin;-M*ones(n,1);DGf_std_min;-M*ones(m,1);zeros(n,1);zeros(2*m,1)];
model.ub  = [vmax;M*ones(n,1);DGf_std_max;M*ones(m,1);ones(n,1);ones(2*m,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(idxNotExch,:);

% Define problem matrix
model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m+2*m);...                      % Sflux*v = 0
    zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n),zeros(n,2*m);...         % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
    zeros(n,nflux),-eye(n),Sthermo',RT*Sthermo',zeros(n),zeros(n,2*m);...          % -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) <= tol
    -Vblock,zeros(n,n+2*m),K*eye(n),zeros(n,2*m);...                               % -v + K*e <= K
    Vblock,zeros(n,n+2*m),-K*eye(n),zeros(n,2*m);...                               % v - K*e <= 0
    zeros(n,nflux),eye(n),zeros(n,2*m),K*eye(n),zeros(n,2*m);...                   % DGr + K*e <= K - delta
    zeros(n,nflux),-eye(n),zeros(n,2*m),-K*eye(n),zeros(n,2*m);...                 % -DGr - K*e <= -delta
    zeros(m,nflux),zeros(m,n),zeros(m,m),eye(m),zeros(m,n),-K*eye(m),zeros(m);...  % ln(x) - K*ub <= ln(xmax)
    zeros(m,nflux),zeros(m,n),zeros(m,m),-eye(m),zeros(m,n),zeros(m),-K*eye(m)]);  % -ln(x) - K*lb <= -ln(xmin)

if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.obj = [zeros(size(model.A,2)-2*m,1); ones(2*m,1)];

% Constraints sense and rhs
model.rhs = [zeros(size(Sflux,1),1);tol*ones(n,1);tol*ones(n,1);K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1);log(xmax);-log(xmin)];
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

% Define optimization parameters
params.outputflag    = 0;
params.IntFeasTol    = 1e-9;
params.MIPGap        = 1e-6;
params.MIPGapAbs     = 1e-12;
params.OptimalityTol = 1e-9;

% Check the feasibility of the problem
sol = gurobi(model,params);

rowList = find(sol.x(end-2*m+1:end)==1);
boundaryList = cell(numel(rowList),1);

for entry=1:numel(rowList)
    row = rowList(entry);
    
    if row > m
        row = row-m;
        rowList(entry) = row;
        boundaryList{entry} = 'lower';
    else
        boundaryList{entry} = 'upper';
    end
    
end

metList = ensemble.mets(rowList);


end

