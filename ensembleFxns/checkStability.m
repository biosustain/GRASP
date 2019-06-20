function maxEigenvalue = checkStability(ensemble,models,strucIdx)
%---------------- Pedro Saa UQ 2018, Marta Matos 2019 ---------------------

% Optimization & simulation parameters
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = numel(ensemble.freeVars);
Sred         = ensemble.Sred;
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsActive);
ix_enz       = ix_mets(end)+1:freeVars;

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
simFlux = feval(kineticFxn,xstep,models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);

% Compute elasticiy matrices
E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))'; % equivalent to imag(simFlux)./1.0e-10 ? 

% Compute Jacobian eigenvalues
%C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
jacobian   = Sred*E_x_abs;
eigenvalues = eig(jacobian);

% Look for positive real eigenvalues
maxEigenvalue = max(real(eigenvalues));

end





