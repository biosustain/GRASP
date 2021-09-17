function findIssuesWithTMFA(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,gurobiModel,options,params)
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
% Problematic metabolites and/or reactions can be missed due to a 
% large K value, especially when using intlinprog (same issue as in
% computeGibbsFreeEnergyRanges). Thus the loop for decreasing K values 
% when the lists of problematic reactions and metabolites are empty.
%
%
% USAGE:
%
%    findIssuesWithTMFA(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,gurobiModel,options,params)
%
% INPUT:
%    ensemble (struct):             model ensemble, see buildEnsemble for fields description
%    DGr_std_min (double vector):   minimum values for standard reaction Gibbs energies
%    DGr_std_max (double vector):   maximum values for standard reaction Gibbs energies
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
%    gurobiModel (struct):          MILP problem structure for gurobi
%    options (struct):              configurations for intlinprog
%    params (struct):               configurations for gurobi
%
%
% .. Authors:
%       - Marta Matos	2020 original code


[metsList, metsBoundaryList] = findProblematicMetabolites(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params);
[rxnsList, rxnsBoundaryList] = findProblematicReactions(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params);

while isempty(metsList) && isempty(rxnsList) && K >= max(abs(vmax), [], 'all')
    n = log10(K);
    K = 10^(n-1);
    [metsList, metsBoundaryList] = findProblematicMetabolites(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params);
    [rxnsList, rxnsBoundaryList] = findProblematicReactions(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params);

end

if isempty(metsList) && isempty(rxnsList) && strcmp(ensemble.LPSolver, 'gurobi')
    gurobi_write(gurobiModel, 'TMFA_problem.lp');
    
    solIIS = gurobi_iis(gurobiModel);
    constraintsIssues = find(solIIS.Arows == 1);
    lbIssues = find(solIIS.lb == 1);
    ubIssues = find(solIIS.ub == 1);

    error(strcat('The TMFA is infeasible. The MILP model has been written in the examples folder so that you can inspect it (TMFA_problem.lp). There might be issues with constraints number: ', num2str(constraintsIssues), ' variable lower bounds number ', num2str(lbIssues), ' or upper bounds number ', num2str(ubIssues), '. '));
else
    if ~isempty(metsList)
        disp(['There is a problem with the metabolite bounds ', strjoin(metsBoundaryList, ', ') ,' for metabolites ', strjoin(metsList, ', '), ' in thermoMets']);
    end
   if ~isempty(rxnsList)
        disp(['There is a problem with the reaction bounds ', strjoin(rxnsBoundaryList, ', '), ' for reactions ', strjoin(rxnsList, ', '), ' in thermoRxns']);
    end
    error(strcat('The TMFA problem is infeasible, see message above. Verify that the standard Gibbs free energy and metabolite concentration values are valid/correct. Note that the problem can also be with the reaction fluxes. Bottom line, for each reaction, fluxes and Gibbs energies need to agree.'));
end
    
end


function [rxnsList, boundaryList] = findProblematicReactions(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params)
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
%    [rxnsList, boundaryList] = findProblematicReactions(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params)
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
model.lb  = [vmin;-M*ones(n,1);DGr_std_min;log(xmin);zeros(3*n,1)];
model.ub  = [vmax;M*ones(n,1);DGr_std_max;log(xmax);ones(3*n,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(ensemble.idxNotExch,:);
Nint = null(Sthermo,'r');
nRowsNint = size(Nint', 1);

% Define problem matrix

model.Aeq = sparse([Sflux,zeros(size(Sflux,1),5*n+m)]);
model.A   = sparse([zeros(n,nflux),eye(n),-eye(n),-RT*Sthermo',zeros(n),-K*eye(n),zeros(n,n);...   % DGr - DGr_std - RT*Sthermo'*ln(x) - K*ubound <= tol
                    zeros(n,nflux),-eye(n),eye(n),RT*Sthermo',zeros(n),zeros(n,n),-K*eye(n);...    % -DGr + DGr_std + RT*Sthermo'*ln(x) - K*lbound <= tol
                    -Vblock,zeros(n,2*n+m),K*eye(n),zeros(n,2*n);...                               % -v + K*e <= K
                    Vblock,zeros(n,2*n+m),-K*eye(n),zeros(n,2*n);...                               % v - K*e <= 0
                    zeros(n,nflux),eye(n),zeros(n,n+m),K*eye(n),zeros(n,2*n);...                   % DGr + K*e <= K - delta
                    zeros(n,nflux),-eye(n),zeros(n,n+m),-K*eye(n),zeros(n,2*n);...                 % -DGr - K*e <= -delta
                    zeros(nRowsNint,nflux),Nint',zeros(nRowsNint,4*n+m);...                        % Nint'*DGr <= tol
                    zeros(nRowsNint,nflux),-Nint',zeros(nRowsNint,4*n+m)]);                        % -Nint'*DGr <= tol

if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.f = [zeros(size(model.A,2)-2*n,1); ones(2*n,1)];

% Constraints sense and rhs
model.beq = zeros(size(Sflux,1),1);
model.b   = [tol*ones(n,1);tol*ones(n,1);K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1);tol*ones(2*nRowsNint,1)];

if ~isempty(ineqConstraints)
    model.b = [model.rhs;ineqConstraints(:,end)];
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
    gurobiModel.vtype(1:(nflux+2*n+m)) = 'C';
    gurobiModel.vtype((nflux+2*n+m+1):end) = 'B';

    % Check the feasibility of the problem
    sol = gurobi(gurobiModel,params);
    if ~strcmp(sol.status,'INFEASIBLE')
        solX = sol.x;
    end


elseif strcmp(ensemble.LPSolver, 'linprog')
    
    model.intcon = [(nflux+2*n+m+1):numel(model.f)];
    [solX,fval] = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);

end


    
        
if (strcmp(ensemble.LPSolver, 'gurobi') && strcmp(sol.status,'INFEASIBLE')) || ...
   (strcmp(ensemble.LPSolver, 'linprog') && isempty(fval))
        
    rxnsList = [];
    boundaryList = [];
    disp('Infeasible solution when trying to find problematic reactions.');    
    return;
end

rowList = find(solX(end-2*n+1:end)==1);
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



function [metsList, boundaryList] = findProblematicMetabolites(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params)
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
%    [metsList, boundaryList] = findProblematicMetabolites(ensemble,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol,options,params)
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
model.lb  = [vmin;-M*ones(n,1);DGr_std_min;-M*ones(m,1);zeros(n,1);zeros(2*m,1)];
model.ub  = [vmax;M*ones(n,1);DGr_std_max;M*ones(m,1);ones(n,1);ones(2*m,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(ensemble.idxNotExch,:);
Nint = null(Sthermo,'r');
nRowsNint = size(Nint', 1);

% Define problem matrix
model.Aeq = sparse([Sflux,zeros(size(Sflux,1),3*n+3*m)]);
model.A  = sparse([zeros(n,nflux),eye(n),-eye(n),-RT*Sthermo',zeros(n),zeros(n,2*m);...            % DGr - DGr_std - RT*Sthermo'*ln(x) <= tol
                   zeros(n,nflux),-eye(n),eye(n),RT*Sthermo',zeros(n),zeros(n,2*m);...             % -DGr + DGr_std + RT*Sthermo'*ln(x) <= tol
                   -Vblock,zeros(n,2*n+m),K*eye(n),zeros(n,2*m);...                                % -v + K*e <= K
                   Vblock,zeros(n,2*n+m),-K*eye(n),zeros(n,2*m);...                                % v - K*e <= 0
                   zeros(n,nflux),eye(n),zeros(n,n+m),K*eye(n),zeros(n,2*m);...                    % DGr + K*e <= K - delta
                   zeros(n,nflux),-eye(n),zeros(n,n+m),-K*eye(n),zeros(n,2*m);...                  % -DGr - K*e <= -delta
                   zeros(m,nflux),zeros(m,2*n),eye(m),zeros(m,n),-K*eye(m),zeros(m);...            % ln(x) - K*ubound <= ln(xmax)
                   zeros(m,nflux),zeros(m,2*n),-eye(m),zeros(m,n),zeros(m),-K*eye(m);...           % -ln(x) - K*lbound <= -ln(xmin)
                   zeros(nRowsNint,nflux),Nint',zeros(nRowsNint,2*n+3*m);...                       % Nint'*DGr <= tol
                   zeros(nRowsNint,nflux),-Nint',zeros(nRowsNint,2*n+3*m)]);                       % -Nint'*DGr <= tol

if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.f = [zeros(size(model.A,2)-2*m,1); ones(2*m,1)];

% Constraints sense and rhs
model.beq = zeros(size(Sflux,1),1);
model.b   = [tol*ones(n,1);tol*ones(n,1);K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1);log(xmax);-log(xmin);tol*ones(2*nRowsNint,1)];
if ~isempty(ineqConstraints)
    model.rhs = [model.rhs;ineqConstraints(:,end)];
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
    gurobiModel.vtype(1:(nflux+2*n+m)) = 'C';
    gurobiModel.vtype((nflux+2*n+m+1):end) = 'B';

    % Check the feasibility of the problem
    sol = gurobi(gurobiModel,params);
    if ~strcmp(sol.status,'INFEASIBLE')
        solX = sol.x;
    end

elseif strcmp(ensemble.LPSolver, 'linprog')

    model.intcon = [(nflux+2*n+m+1):numel(model.f)];

    [solX,fval] = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);

end


if (strcmp(ensemble.LPSolver, 'gurobi') && strcmp(sol.status,'INFEASIBLE')) || ...
   (strcmp(ensemble.LPSolver, 'linprog') && isempty(fval))
   
    metsList = [];
    boundaryList = [];
    disp('Infeasible solution when trying to find problematic metabolites.');
    return;
end
    
rowList = find(solX(end-2*m+1:end)==1);
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