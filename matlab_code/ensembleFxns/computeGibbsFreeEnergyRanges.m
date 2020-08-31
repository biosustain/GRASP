function [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints)
% Applies Thermodynamic Flux Balance Analysis (TMFA) to find
% thermodynamically feasible ranges for reactions fluxes, metabolite
% concentrations, and reactions Gibbs energies.
%
% It is based on https://arxiv.org/pdf/1803.04999.pdf
%
% The standard Gibbs energies are converted into metabolite
% formation energies, so that the dG values of reactions 
% involved in a loop sum to zero.
% Note that those formation energies are compatible with the
% given standard Gibbs energies, but are not necessarily
% similar to formation energies obtained from e.g. eQuilibrator.
%
%
% USAGE:
%
%    [v_range,DGr_range,DGf_std_range,lnx_range,x0] = computeGibbsFreeEnergyRanges(ensemble,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,ineqConstraints)
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
%
% OUTPUT:
%    v_range (double matrix):        range of feasible reaction fluxes
%    DGr_range (double matrix):	     range of feasible Gibbs energies
%    DGf_std_range (double matrix):	 range of feasible formation energies
%    lnx_range (double matrix):      range of feasible metabolite concentrations (ln)
%    x0 (double matrix):             initial point to be used for the sampling with hit-and-run method
%
% .. Authors:
%       - Pedro Saa                 2016 original code
%       - Marta Matos	            2018 added error messanges and findProblematicReactions
%       - Pedro Saa and Marta Matos	2020  modified the TMFA MILP to use formation energies

Sthermo   = ensemble.Sthermo;
Sflux     = ensemble.Sflux;

DGf_std_max = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_max;
DGf_std_min = pinv(Sthermo*Sthermo')*Sthermo*DGr_std_min;
for ix = 1:numel(DGf_std_max)
    if (DGf_std_min(ix)>DGf_std_max(ix))
        DGf_std_temp    = DGf_std_max(ix);
        DGf_std_max(ix) = DGf_std_min(ix);
        DGf_std_min(ix) = DGf_std_temp;
    end
end


% Build the adapted TMFA problem
K       = 1e12;
delta   = 1e-6;
RT      = 8.314*298.15/1e3;  % [kJ/mol]
M       = 1e5;
tol     = 1e-10;

% Define bounds
[m,n]     = size(Sthermo);
[~,nflux] = size(Sflux);
model.lb  = [vmin;-M*ones(n,1);DGf_std_min;log(xmin);zeros(n,1)];
model.ub  = [vmax;M*ones(n,1);DGf_std_max;log(xmax);ones(n,1)];
Vblock    = eye(nflux);
Vblock    = Vblock(ensemble.idxNotExch,:);

% Define problem matrix
model.A  = sparse([Sflux,zeros(size(Sflux,1),2*n+2*m);...         % Sflux*v = 0
    zeros(n,nflux),eye(n),-Sthermo',-RT*Sthermo',zeros(n);...     % DGr - Sthermo'*DGf_std - RT*Sthermo'*ln(x) <= tol
    zeros(n,nflux),-eye(n),Sthermo',RT*Sthermo',zeros(n);...      % -DGr + Sthermo'*DGf_std + RT*Sthermo'*ln(x) <= tol
    -Vblock,zeros(n,n+2*m),K*eye(n);...                           % -v + K*e <= K
    Vblock,zeros(n,n+2*m),-K*eye(n);...                           % v - K*e <= 0
    zeros(n,nflux),eye(n),zeros(n,2*m),K*eye(n);...               % DGr + K*e <= K - delta
    zeros(n,nflux),-eye(n),zeros(n,2*m),-K*eye(n)]);              % -DGr - K*e <= -delta

if ~isempty(ineqConstraints)
    p = size(ineqConstraints,1);
    model.A = sparse([model.A;zeros(p,nflux+n+m),ineqConstraints(:,1:end-1),zeros(p,n)]);
end

% Objective function
model.obj = zeros(size(model.A,2),1);

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

if strcmp(sol.status,'INFEASIBLE')
    disp('TMFA problem is infeasible.')
    findIssuesWithTMFA(ensemble,model,DGf_std_min,DGf_std_max,vmin,vmax,xmin,xmax,ineqConstraints,K,RT,delta,M,tol);
end


% Run improved TMFA
v_range   = zeros(nflux,2);
DGr_range = zeros(n,2);
lnx_range = zeros(m,2);
DGf_std_range = zeros(m,2);

xprev_min = zeros(nflux+n+2*m,1);
xprev_max = zeros(nflux+n+2*m,1);

for ix = 1:nflux+n+2*m
    model.obj(ix)    = 1;
    
    model.modelsense = 'min';
    solmin = gurobi(model,params);
    if ~strcmp(solmin.status,'OPTIMAL')
        error(strcat('Check your metabolite concentration ranges in thermoMets. In particular for minimum values set to 0 change them to a low number like 10^15.'));
    end
    xprev_min(:,ix) = solmin.x(1:nflux+n+2*m);                            % save previous min solution for later
    
    model.modelsense = 'max';
    solmax = gurobi(model,params);
    xprev_max(:,ix) = solmax.x(1:nflux+n+2*m);                            % save previous max solution for later
    
    if (ix<=nflux)
        v_range(ix,:) = [solmin.objval,solmax.objval];                    % [umol or mmol/gdcw/h]
    elseif (ix<=nflux+n)
        DGr_range(ix-nflux,:) = [solmin.objval,solmax.objval];            % [kJ/mol]
    elseif (ix<=nflux+n+m)
        DGf_std_range(ix-nflux-n,:) = [solmin.objval,solmax.objval];      % [kJ/mol]
    else
        lnx_range(ix-nflux-n-m,:)   = [solmin.objval,solmax.objval];      % log([M])
    end
    model.obj(ix) = 0;
end


% Return initial point
x0 = mean([xprev_min,xprev_max],2);

Nint = null(Sthermo,'r');
loopCondition = Nint' * x0(nflux+1:nflux+n);

% Check that the initial point is feasible
Acheck = full(model.A(1:(size(Sflux,1)+2*n),1:(nflux+n+2*m)));              % Generate checking matrix
if all(abs(Acheck*x0)<1e-6)
    if ~isempty(Nint)
        if max(abs(loopCondition)) <= 1e-6
            disp('The initial point obtained from TMFA is feasible and valid for starting the sampler.');
        else
            error('The initial point obtained from TMFA is not thermodynamically feasible.')
        end
    else
        disp('The initial point obtained from TMFA is feasible and valid for starting the sampler.');
    end            

    if all(sign(xprev_min(1:nflux))&sign(xprev_max(1:nflux)))    
        disp('Flux directions are consistent.');
    else 
        error('The flux directions are not consistent. Please make sure that both the lower and upper bound of the flux ranges (fluxMean - 2*fluxStd and fluxMean + 2*fluxStd, respectively) are either positive or negative.');
    end
    
else
    error('The initial point obtained from TMFA is not thermodynamically feasible.')
end

end