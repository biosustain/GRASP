function [rowList, dGList] = findProblematicReactions(gurobiModel,params,model,options,DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,rxnNames, solver)
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


DGr_std_min_orig = DGr_std_min;
DGr_std_max_orig = DGr_std_max;

DGr_std_min = ones(size(DGr_std_min)) * -100;
DGr_std_max = ones(size(DGr_std_max)) * 100;

row = 1;
rowList = [];
dGList = [];

sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,gurobiModel,params,model,options,solver);

while ((strcmp(solver, 'gurobi') && strcmp(sol.status,'OPTIMAL'))  || (strcmp(solver, 'linprog') && ~isempty(sol))) ... 
        && row <= numel(DGr_std_min)
    
    DGr_std_min(row) = DGr_std_min_orig(row);
    DGr_std_max(row) = DGr_std_max_orig(row);
    
    sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,gurobiModel,params,model,options,solver);

    if (strcmp(solver, 'gurobi') && strcmp(sol.status,'INFEASIBLE')) || (strcmp(solver, 'linprog') && isempty(sol))
        DGr_std_min(row) = -100;
        DGr_std_max(row) = 100;
        rowList = [rowList, rxnNames(row)];
        dGList = [dGList, strcat('[', num2str(DGr_std_min_orig(row)),',', num2str(DGr_std_max_orig(row)),']')];
        
        sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,gurobiModel,params,model,options,solver);

    end 
    row = row +1;
        
end

end


function [sol] = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,gurobiModel,params,model,options,solver)


if strcmp(solver, 'gurobi')
    gurobiModel.rhs = [zeros(size(Sflux,1),1);DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];

    if ~isempty(ineqConstraints)
        gurobiModel.rhs = [gurobiModel.rhs;ineqConstraints(:,end)];
    end

    sol = gurobi(gurobiModel,params);

elseif strcmp(solver, 'linprog')
    model.b = [DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];
    
    if ~isempty(ineqConstraints)
        model.b = [model.b;ineqConstraints(:,end)];
    end
    
    [x,sol] = intlinprog(model.f, model.intcon, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options);
    
end
end

