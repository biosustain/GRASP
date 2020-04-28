function isModelValid = checkStability(ensemble,models,strucIdx,eigThreshold)
% Given a model, checks if the largest real part of the jacobian
% eigenvalues is greater than the given threshold (*eigThreshold*), if so 
% the model is considered invalid.
%
%
% USAGE:
%
%    isModelValid = checkStability(ensemble, models, strucIdx, eigThreshold)
%
% INPUT:
%    ensemble (struct):       model ensemble, see buildEnsemble for fields description
%    models (struct):         model, see initialSampler for fields description
%    strucIdx (int):          ID of the model structure considered
%    eigThreshold (double):	  threshold for positive eigenvalues' real part
%
% OUTPUT:
%    isModelValid (logical):	whether or not the model is valid
%
% .. Authors:
%       - Pedro Saa     2018 original code
%       - Marta Matos	2019 refactored code to be used in initialSampler

% Optimization & simulation parameters
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = numel(ensemble.freeVars);
Sred         = ensemble.Sred;
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsActive);
ix_enz       = ix_mets(end)+1:freeVars;
xconst       = ones(numel(ensemble.metsFixed), numel(ix_mets));

% Main loop
hstep = 1e-10;              % Step size for control coefficient computations
xopt = ones(freeVars,1);

% Define reference state
xref = xopt(ix_mets);
Eref = xopt(ix_enz);

% Define step length to perturb metabolite concentrations
hstep_x = hstep*xref;
xmets   = repmat(xref,1,numel(xref)) + 1i*diag(hstep_x);
xenz    = repmat(Eref,1,numel(xref));
xstep   = [xmets;xenz];

% Simulate flux for metabolite perturbation
%simFlux = feval(kineticFxn,xstep,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0);
simFlux = feval(kineticFxn,xstep,xconst,models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);

% Compute elasticiy matrices
E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))'; % equivalent to imag(simFlux)./1.0e-10 ? 

% Compute Jacobian eigenvalues
jacobian   = Sred*E_x_abs;
eigenvalues = eig(jacobian);

% Look for positive real eigenvalues
maxEigenvalue = max(real(eigenvalues));

isModelValid = true;

if maxEigenvalue > eigThreshold
    isModelValid = false;
end

end





