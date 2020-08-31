function findIssuesWithTMFA(ensemble,model,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
% Finds out why the TMFA problem is infeasible.
%
% This is done by first trying to find which reactions' standard 
% Gibbs energies or metabolite concentrations are causing the 
% problem to be infeasible.
% If no reactions nor metabolites are deemed to be the cause and
% if the LP solver is set to gurobi, the MILP is written into 
% a text file and the function `gurobi_iis` is used to find
% which constraints/variable bounds are causing the issue.
%
%
% USAGE:
%
%    findIssuesWithTMFA(ensemble,model,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
%
% INPUT:
%    ensemble (struct):             model ensemble, see buildEnsemble for fields description
%    model (struct):                MILP problem structure
%    DGf_std_min (double vector):   minimum values for metabolite formation energies
%    DGf_std_max (double vector):   maximum values for metabolite formation energies
%    vmin (double vector):          minimum values for reaction fluxes
%    vmax (double vector):          maximum values for reaction fluxes
%    xmin (double vector):	        minimum values for metabolite concentrations
%    xmax (double vector):          maximu values for metabolite concentrations
%    ineqConstraints (double):	    inequality constraints
%    K (double):                    large value used in the MILP problem
%    RT (double):                   value for the gas constant multiplied by the temperature (in Kelvin)
%    delta (double):                low value to use in the MILP problem
%    M (double):                    large value to use in the MILP problem
%    tol (double):                  tolerance value to use when comparing results to zero.
%
%
% .. Authors:
%       - Marta Matos	2020 original code


[metsList, metsBoundaryList] = findProblematicMetabolites(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol);

[rxnsList, rxnsBoundaryList] = findProblematicReactions(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol);

if isempty(metsList) && isempty(rxnsList)
    
    if strcmp(ensemble.LPSolver, 'gurobi')
        gurobi_write(model, 'TMFA_problem.lp');
        solIIS = gurobi_iis(model);
        constraintsIssues = find(solIIS.Arows == 1);
        lbIssues = find(solIIS.lb == 1);
        ubIssues = find(solIIS.ub == 1);
        
        error(strcat('The TMFA is infeasible. The MILP model has been written in the examples folder so that you can inspect it (TMFA_problem.lp). There might be issues with constraints number: ', num2str(constraintsIssues), ' variable lower bounds number ', num2str(lbIssues), ' or upper bounds number ', num2str(ubIssues), '. '));
    end
else
    error(strcat('The TMFA problem is infeasible. This can be a problem with the metabolite bounds ', strjoin(metsBoundaryList, ', ') ,' for metabolites ', strjoin(metsList, ', '), ' in thermoMets, or a problem with the reaction bounds ', strjoin(rxnsBoundaryList, ', '), ' for reactions ', strjoin(rxnsList, ', '), ' in thermoRxns.Verify that the standard Gibbs free energy and metabolite concentration values are valid/correct. Note that the problem can also be with the reaction fluxes. Bottom line, for each reaction, fluxes and Gibbs energies need to agree.'));
end
    
end


function [rxnsList, boundaryList] = findProblematicReactions(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
% Finds the reactions whose standard Gibbs energy bounds make the 
% TMFA problem infeasible.
%
% This is done by adding a binary variable to the TMFA problem while 
% minimizing their sum. That way, if the binary variable is one 
% instead of zero, the respective bound needs to be relaxed.
%
% The constraints of the modified TMFA problem are the following:
%
% Sflux*v = 0
% DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) + K*ubound <= tol
% -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) - K*lbound <= tol
% -v + K*e <= K
% v - K*e <= 0
% DGr + K*e <= K - delta
% -DGr - K*e <= -delta
%
% where lbound and ubound are the binary variables.
%
%
% USAGE:
%
%    [rxnsList, boundaryList] = findProblematicReactions(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
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
%    K (double):                    large value used in the MILP problem
%    RT (double):                   value for the gas constant multiplied by the temperature (in Kelvin)
%    delta (double):                low value to use in the MILP problem
%    M (double):                    large value to use in the MILP problem
%    tol (double):                  tolerance value to use when comparing results to zero.
%
% OUTPUT:
%    rxnsList (int vector):	        names of reactions that make the TMFA problem infeasible
%    boundaryList (double vector):  dG values of reactions that make the TMFA problem infeasible
%
% .. Authors:
%       - Marta Matos	2019 original code


% Define bounds
Sthermo   = ensemble.Sthermo;
Sflux     = ensemble.Sflux;
[m,n]     = size(Sthermo);
[~,nflux] = size(Sflux);
model.lb  = [vmin;-M*ones(n,1);DGf_std_min;log(xmin);zeros(3*n,1)];
model.ub  = [vmax;M*ones(n,1);DGf_std_max;log(xmax);ones(3*n,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(ensemble.idxNotExch,:);

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

if ~strcmp(sol.status,'OPTIMAL')
    rxnsList = [];
    boundaryList = [];
    disp('Infeasible solution when trying to find problematic reactions.');    
    return;
end

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



function [metsList, boundaryList] = findProblematicMetabolites(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
% Finds the metabolites whose concentration bounds make the 
% TMFA problem infeasible.
%
% This is done by adding a binary variable to the TMFA problem while 
% minimizing their sum. That way, if the binary variable is one 
% instead of zero, the respective bound needs to be relaxed.
%
% The constraints of the modified TMFA problem are the following:
%
% Sflux*v = 0
% DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
% -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) <= tol
% -v + K*e <= K
% v - K*e <= 0
% DGr + K*e <= K - delta
% -DGr - K*e <= -delta
% ln(x) - K*ubound <= ln(xmax)
% -ln(x) - K*lbound <= -ln(xmin)
%
% where lbound and ubound are the binary variables.
%
%
% USAGE:
%
%    [metsList, boundaryList] = findProblematicMetabolites(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol)
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
%    K (double):                    large value used in the MILP problem
%    RT (double):                   value for the gas constant multiplied by the temperature (in Kelvin)
%    delta (double):                low value to use in the MILP problem
%    M (double):                    large value to use in the MILP problem
%    tol (double):                  tolerance value to use when comparing results to zero.
%
% OUTPUT:
%    metsList (int vector):	        names of metabolites that make the LP infeasible
%    boundaryList (double vector):  dG values of reactions that make the LP infeasible
%
% .. Authors:
%       - Marta Matos	2020 original code



% Define bounds
Sthermo   = ensemble.Sthermo;
Sflux     = ensemble.Sflux;
[m,n]     = size(Sthermo);
[~,nflux] = size(Sflux);
model.lb  = [vmin;-M*ones(n,1);DGf_std_min;-M*ones(m,1);zeros(n,1);zeros(2*m,1)];
model.ub  = [vmax;M*ones(n,1);DGf_std_max;M*ones(m,1);ones(n,1);ones(2*m,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(ensemble.idxNotExch,:);

% Define problem matrix
model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m+2*m);...                      % Sflux*v = 0
    zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n),zeros(n,2*m);...         % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
    zeros(n,nflux),-eye(n),Sthermo',RT*Sthermo',zeros(n),zeros(n,2*m);...          % -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) <= tol
    -Vblock,zeros(n,n+2*m),K*eye(n),zeros(n,2*m);...                               % -v + K*e <= K
    Vblock,zeros(n,n+2*m),-K*eye(n),zeros(n,2*m);...                               % v - K*e <= 0
    zeros(n,nflux),eye(n),zeros(n,2*m),K*eye(n),zeros(n,2*m);...                   % DGr + K*e <= K - delta
    zeros(n,nflux),-eye(n),zeros(n,2*m),-K*eye(n),zeros(n,2*m);...                 % -DGr - K*e <= -delta
    zeros(m,nflux),zeros(m,n),zeros(m,m),eye(m),zeros(m,n),-K*eye(m),zeros(m);...  % ln(x) - K*ubound <= ln(xmax)
    zeros(m,nflux),zeros(m,n),zeros(m,m),-eye(m),zeros(m,n),zeros(m),-K*eye(m)]);  % -ln(x) - K*lbound <= -ln(xmin)

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

if ~strcmp(sol.status,'OPTIMAL')
    metsList = [];
    boundaryList = [];
    disp('Infeasible solution when trying to find problematic metabolites.');
    return;
end
    
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

metsList = ensemble.mets(rowList);


end