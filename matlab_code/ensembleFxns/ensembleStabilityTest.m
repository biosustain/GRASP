function stabilityRes = ensembleStabilityTest(ensemble,eigThreshold,strucIdx)
%
% Takes in a model ensemble, calculates the jacobian and respective
% eigenvalues for each model.
% If any eigenvalue's real part is higher than eighThreshold, the model
% is considered unstable.
%
%---------------- Pedro Saa UQ 2018, Marta Matos 2019 ---------------------

if nargin<3
    strucIdx = 1;
    if ensemble.populations(end).strucIdx(1)==0
        ensemble.populations(end).strucIdx = ones(numel(ensemble.populations(end).strucIdx),1);
    end
end

% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
numModels = size(ensemble.populations.models, 2);
if numModels > numel(particleIdx) 
    numModels   = numel(particleIdx);
end

% Optimization & simulation parameters
fixedExchs   = ensemble.fixedExch;
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = numel(ensemble.freeVars);
Sred         = ensemble.Sred;
kinInactRxns = ensemble.kinInactRxns;
subunits     = ensemble.subunits{strucIdx};
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsActive);
ix_enz       = ix_mets(end)+1:freeVars;
metNames     = ensemble.mets(ensemble.metsActive);
rxnNames     = ensemble.rxns;

% Check sampler mode to determine the numer of conditions
if ~strcmpi(ensemble.sampler,'ORACLE')
    nCondition   = size(ensemble.expFluxes,2)+1;
else
    nCondition = 1;
end

% Main loop
hstep = 1e-10;              % Step size for control coefficient computations

stabilityRes.unstableModels = [];
for ix = 1:nCondition
   
    for jx = 1:numModels
        stabilityRes.posEig{jx} = [];
        %mcaResults.enzNames = rxnNames;
        model = ensemble.populations(end).models(particleIdx(jx));
        if ix == 1
            xopt = ones(freeVars,1);
            vref = feval(kineticFxn,xopt,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0);
        else
            xopt = ensemble.populations(end).xopt{particleIdx(jx)}(:,ix-1);
            vref = ensemble.populations(end).simFluxes{particleIdx(jx)}(:,ix-1);
        end
        
        % Define reference state
        xref = xopt(ix_mets);
        Eref = xopt(ix_enz);
        
        % Define step length to perturb metabolite concentrations
        hstep_x = hstep*xref;
        xmets   = repmat(xref,1,numel(xref)) + 1i*diag(hstep_x);
        xenz    = repmat(Eref,1,numel(xref));
        xstep   = [xmets;xenz];
        
        % Simulate flux for metabolite perturbation
        simFlux = feval(kineticFxn,xstep,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0);
        
        % Compute elasticiy matrices
        E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))'; % equivalent to imag(simFlux)./1.0e-10 ? 
                      
        % Compute Jacobian eigenvalues
        %C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
        jacobian   = Sred*E_x_abs;
        eigenvalues = eig(jacobian);
        
        % Look for positive real eigenvalues
        indList = find(real(eigenvalues) > eigThreshold);
        posEigMets = ensemble.mets(ensemble.metsActive(indList));
                
        if size(indList, 1) > 0
            stabilityRes.posEig{jx} = {posEigMets, real(eigenvalues(indList))};
            stabilityRes.unstableModels = [stabilityRes.unstableModels, jx];
        end
        
    end
    
    disp(['*** Out of ', num2str(numModels), ' models, ', num2str(size(stabilityRes.unstableModels, 2)), ' models are unstable ***']);
end





