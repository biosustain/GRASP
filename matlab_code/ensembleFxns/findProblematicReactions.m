function [rowList, dGList] = findProblematicReactions(model,params,DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,sol,rxnNames)
%--------------------------------------------------------------------------
%
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
% solution to be infeasible when their dG bounds are the original.% 
%
%--------------------- Marta Matos 2019 -----------------------------------


DGr_std_min_orig = DGr_std_min;
DGr_std_max_orig = DGr_std_max;

DGr_std_min = ones(size(DGr_std_min)) * -100;
DGr_std_max = ones(size(DGr_std_max)) * 100;

row = 1;
rowList = [];
dGList = [];

sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params);

while strcmp(sol.status,'OPTIMAL') && row <= numel(DGr_std_min)
    
    DGr_std_min(row) = DGr_std_min_orig(row);
    DGr_std_max(row) = DGr_std_max_orig(row);
    
    sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params);

    if strcmp(sol.status,'INFEASIBLE')
        DGr_std_min(row) = -100;
        DGr_std_max(row) = 100;
        rowList = [rowList, rxnNames(row)];
        dGList = [dGList, strcat('[', num2str(DGr_std_min_orig(row)),',', num2str(DGr_std_max_orig(row)),']')];
        
        sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params);

    end 
    row = row +1;
        
end

end


function [sol] = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params)

model.rhs = [zeros(size(Sflux,1),1);DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];
if ~isempty(ineqConstraints)
    model.rhs = [model.rhs;ineqConstraints(:,end)];
end

sol = gurobi(model,params);

end

