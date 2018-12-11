function [row_list, dg_list] = findProblematicReactions(model,params,DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,sol)

DGr_std_min_orig = DGr_std_min;
DGr_std_max_orig = DGr_std_max;

DGr_std_min = ones(size(DGr_std_min)) * -100;
DGr_std_max = ones(size(DGr_std_max)) * 100;

row = 1;
row_list = [];
dg_list = [];

sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params);

while strcmp(sol.status,'OPTIMAL') && row <= numel(DGr_std_min)
    
    DGr_std_min(row) = DGr_std_min_orig(row);
    DGr_std_max(row) = DGr_std_max_orig(row);
    
    sol = runOptimization(DGr_std_min,DGr_std_max,K,delta,n,Sflux,ineqConstraints,model,params);

    if strcmp(sol.status,'INFEASIBLE')
        DGr_std_min(row) = -100;
        DGr_std_max(row) = 100;
        row_list = [row_list, row];
        dg_list = [dg_list, strcat('[', num2str(DGr_std_min_orig(row)),',', num2str(DGr_std_max_orig(row)),']')];
        
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

