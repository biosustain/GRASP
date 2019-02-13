function mcaResults = controlAnalysis(ensemble,strucIdx)
%---------------- Pedro Saa UQ 2018----------------------------------------

if nargin<2
    strucIdx = 1;
    if ensemble.populations(end).strucIdx(1)==0
        ensemble.populations(end).strucIdx = ones(numel(ensemble.populations(end).strucIdx),1);
    end
end

% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
numModels   = numel(particleIdx);

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
for ix = 1:nCondition
    mcaResults.xControl{ix}    = [];
    mcaResults.xControlAvg{ix} = 0;
    mcaResults.vControl{ix}    = [];
    mcaResults.vControlAvg{ix} = 0;
    mcaResults.xcounter{ix}    = 0;
    mcaResults.vcounter{ix}    = 0;
    mcaResults.E_x_nor{ix}     = [];
    
    for jx = 1:numModels
        mcaResults.enzNames = rxnNames;
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
        
        % Normalize elasticity matrices
        E_x_nor   = diag(vref.^-1)*E_x_abs*diag(xref);        
              
        % Compute control coefficients
        C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
        C_x       = diag(xref.^-1)*C_x_abs*diag(vref);
        C_v       = eye(numel(vref)) + E_x_nor*C_x;

        % Save control coefficients only if the result is accurate
        if all(abs(sum(C_x,2))<1e-5)
            mcaResults.xControl{ix}    = [mcaResults.xControl{ix}; C_x];
            mcaResults.xControlAvg{ix} = mcaResults.xControlAvg{ix} + C_x;
            mcaResults.xcounter{ix}    = mcaResults.xcounter{ix} + 1;
        end
        if all(abs(sum(C_v,2))-1<1e-5)
            mcaResults.vControl{ix}    = [mcaResults.vControl{ix}; C_v];
            mcaResults.vControlAvg{ix} = mcaResults.vControlAvg{ix} + C_v;
            mcaResults.vcounter{ix}    = mcaResults.vcounter{ix} + 1;
            mcaResults.E_x_nor{ix}     = [mcaResults.E_x_nor{ix}; E_x_nor];
        end
    end
    
    % Determine expectancy for control coefficients
    mcaResults.xControlAvg{ix} = mcaResults.xControlAvg{ix}/mcaResults.xcounter{ix};
    mcaResults.vControlAvg{ix} = mcaResults.vControlAvg{ix}/mcaResults.vcounter{ix};
    
end

end
