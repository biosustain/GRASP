function [DGr_rng,xrng,vrng] = TFBA(S,DGr_std,vmin,vmax,xmin,xmax,intMets)
% Thermodynamic-based Flux Balance Analysis
% 
% ------------------ Pedro Saa 2018 ---------------------------------------
%
% Build the adapted TMFA problem (there is one reversible rxn)
K       = 1e4;
delta   = 1e-6;
RT      = 8.314*298.15/1e3;  % [kJ/mol]

% Define bounds
[m,n]    = size(S);
Sint     = S(intMets,:);
model.lb = [vmin;-K*ones(n,1);log(xmin);zeros(n,1)];
model.ub = [vmax;K*ones(n,1);log(xmax);ones(n,1)];

% Define problem matrix
model.A  = sparse([Sint,zeros(size(Sint,1),2*n+m);...       % Sint*v = 0
    zeros(n),eye(n),-RT*S',zeros(n);...                     % DGr - RT*Sint*ln(x) = DGr_std
    -eye(n),zeros(n,n+m),K*eye(n);...                       % -v + K*e <= K
    eye(n),zeros(n,n+m),-K*eye(n);...                       % v - K*e <= 0
    zeros(n),eye(n),zeros(n,m),K*eye(n);...                 % DGr + K*e <= K - delta
    zeros(n),-eye(n),zeros(n,m),-K*eye(n)]);                % -DGr - K*e <= -delta

% Objective function
model.obj = zeros(size(model.A,2),1);

% Constraints sense and rhs
model.rhs   = [zeros(size(Sint,1),1);DGr_std;K*ones(n,1);zeros(n,1);...
    (K-delta)*ones(n,1);-delta*ones(n,1)];
model.sense = blanks(numel(model.rhs));
for ix = 1:numel(model.rhs)
    if (ix<=size(Sint,1)+n)
        model.sense(ix) = '=';
    else
        model.sense(ix) = '<';
    end
end

% Variable type definition
model.vtype = blanks(numel(model.obj));
for ix = 1:numel(model.obj)
    if (ix<=2*n+m)
        model.vtype(ix) = 'C';
    else
        model.vtype(ix) = 'B';
    end
end

% Define optimization parameters
params.outputflag = 0;

% Run improved TMFA
vrng    = zeros(n,2);
DGr_rng = zeros(n,2);
xrng    = zeros(m,2);
for ix = 1:2*n+m
    model.obj(ix)    = 1;
    model.modelsense = 'min';
    solmin           = gurobi(model,params);
    model.modelsense = 'max';
    solmax           = gurobi(model,params);
    if (ix<=n)
        vrng(ix,:)      = [solmin.objval,solmax.objval];          % [mmol/gdcw/h]
    elseif (ix<=2*n)
        DGr_rng(ix-n,:) = [solmin.objval,solmax.objval];          % [kJ/mol]
    else
        xrng(ix-2*n,:)  = [solmin.objval,solmax.objval];          % [M]
    end
    model.obj(ix) = 0;
end
xmax = max(xrng(:));  % robust calculation of exp(log(x))
xrng = xrng - xmax;
xrng = exp(xrng);
xrng = xrng*exp(xmax);
