function [DGr_rng,xrng,vrng] = computeGibbsFreeEnergyRanges(Sflux,Sthermo,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints)
% Thermodynamic-based Flux Balance Analysis
% 
% ------------------ Pedro Saa, Marta Matos 2018 ---------------------------------------
%
if (nargin<10); ineqConstraints = []; end

% Build the adapted TMFA problem
K       = 1e4;
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
model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+m);...           % Sflux*v = 0
    zeros(n,nflux),eye(n),-RT*Sthermo',zeros(n);...               % DGr - RT*Sthermo*ln(x) <= DGr_std_max
    zeros(n,nflux),-eye(n),RT*Sthermo',zeros(n);...               % -DGr + RT*Sthermo*ln(x) <= -DGr_std_min
    -Vblock,zeros(n,n+m),K*eye(n);...                             % -v + K*e <= K
    Vblock,zeros(n,n+m),-K*eye(n);...                             % v - K*e <= 0
    zeros(n,nflux),eye(n),zeros(n,m),K*eye(n);...                 % DGr + K*e <= K - delta
    zeros(n,nflux),-eye(n),zeros(n,m),-K*eye(n)]);                % -DGr - K*e <= -delta
if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.obj = zeros(size(model.A,2),1);

% Constraints sense and rhs
model.rhs = [zeros(size(Sflux,1),1);DGr_std_max;-DGr_std_min;K*ones(n,1);zeros(n,1);(K-delta)*ones(n,1);-delta*ones(n,1)];
if ~isempty(ineqConstraints)
    model.rhs = [model.rhs;ineqConstraints(:,end)];
end
model.sense = blanks(numel(model.rhs));
for ix = 1:numel(model.rhs)
    if (ix<=size(Sflux,1))
        model.sense(ix) = '=';
    else
        model.sense(ix) = '<';
    end
end

% Variable type definition
model.vtype = blanks(numel(model.obj));
for ix = 1:numel(model.obj)
    if (ix<=nflux+n+m)
        model.vtype(ix) = 'C';
    else
        model.vtype(ix) = 'B';
    end
end

% Define optimization parameters
params.outputflag = 0;

% Check the feasibility of the problem
sol = gurobi(model,params);
assert(~strcmp(sol.status,'INFEASIBLE'),...
    'The TMFA problem is infeasible. Verify that the standard Gibbs free energy and metabolite concentration values are valid/correct.');

% Run improved TMFA
vrng    = zeros(nflux,2);
DGr_rng = zeros(n,2);
xrng    = zeros(m,2);
for ix = 1:nflux+n+m
    model.obj(ix)    = 1;
    model.modelsense = 'min';
    solmin           = gurobi(model,params);
    model.modelsense = 'max';
    solmax           = gurobi(model,params);
    if (ix<=nflux)
        vrng(ix,:)   = [solmin.objval,solmax.objval];             % [umol or mmol /gdcw/h]
    elseif (ix<=nflux+n)
        DGr_rng(ix-nflux,:) = [solmin.objval,solmax.objval];      % [kJ/mol]
    else
        xrng(ix-nflux-n,:)  = [solmin.objval,solmax.objval];      % exp([M])
    end
    model.obj(ix) = 0;
end
xmax = max(xrng(:));  % robust calculation of exp(log(x))
xrng = xrng - xmax;
xrng = exp(xrng);
xrng = xrng*exp(xmax);