function mcaResults = controlAnalysis(ensemble,saveResMatrices,conditionI,strucIdx)
% Do Metabolic Control Analysis for all models in the given ensemble and
% return the average flux and concentration control coefficients across the 
% ensemble.
%
% It can also return the control coefficient and elasticity matrices for
% each model in the ensemble if *saveResMatrices* is set to true. But keep 
% in mind that for large models this can use a lot of RAM and CPU.
%
% The implementation is based on 
%
% - Computational Models of Metabolism  Stability and Regulation in 
%   Metabolic Networks, section 1.7.2 on Metabolic Control Analysis
%   (https://doi.org/10.1002/9780470475935.ch3).
%
%
% USAGE:
%
%    mcaResults = controlAnalysis(ensemble, saveResMatrices, strucIdx)
%
% INPUT:
%    ensemble (struct):           model ensemble, see buildEnsemble for fields description
%    saveResMatrices (logical):   whether or not to save the elasticity and control coefficient matrices for all models
%
% OPTIONAL INPUT:
%    conditionI (int):  condition for which the user wants for run MCA when having multiple coniditions, by
%                       default MCA is run for all conditions
%    strucIdx (int):	index of the model structure considered
%
% OUTPUT:
%    mcaResults (struct):	MCA results
%
%               * xControlAvg (*cell*)   : average concentration control coefficient for the model ensemble
%               * vControlAvg (*cell*)   : average flux control coefficient for the model ensemble
%               * xcounter (*cell*)      : number of models in the average concentration control coefficient calculation
%               * vcounter (*cell*)      : number of models in the average flux control coefficient calculation
%               * xControl (*cell*)      : concentration control coefficient matrix for each model
%               * vControl (*cell*)      : flux control coefficient matrix for each model
%               * E_x_nor (*cell*)       : normalized elasticity matrix for each model
%
% .. Authors:
%       - Pedro Saa     2018 original code
%       - Marta Matos   2019 refactored code

if nargin<4
    strucIdx = 1;
    if ensemble.populations(end).strucIdx(1)==0
        ensemble.populations(end).strucIdx = ones(numel(ensemble.populations(end).strucIdx),1);
    end
end

% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
numModels   = numel(ensemble.populations.models);

% Optimization & simulation parameters
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = numel(ensemble.freeVars);
Sred         = ensemble.Sred;
kinInactRxns = ensemble.kinInactRxns;
subunits     = ensemble.subunits{strucIdx};
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsBalanced);
ix_enz       = ix_mets(end)+1:freeVars;
metNames     = ensemble.mets(ensemble.metsBalanced);
rxnNames     = ensemble.rxns;

% Check sampler mode to determine the numer of conditions
if ~strcmpi(ensemble.sampler,'GRASP')
    nCondition   = size(ensemble.expFluxes,2)+1;
else
    nCondition = 1;
end

if nargin<3
    startCondition = 1;
    endCondition = nCondition;
else
    startCondition = conditionI;
    endCondition = conditionI;
end

% Main loop
hstep = 1e-10;              % Step size for control coefficient computations
for ix = startCondition:endCondition
    if saveResMatrices
        mcaResults.xControl{ix}    = [];
        mcaResults.vControl{ix}    = [];
        mcaResults.E_x_nor{ix}     = [];
    end
    
    mcaResults.xControlAvg{ix} = 0;
    mcaResults.vControlAvg{ix} = 0;
    mcaResults.xcounter{ix}    = 0;
    mcaResults.vcounter{ix}    = 0;
    
    for jx = 1:numModels
        disp(['Model: ', num2str(jx)]);
        
        if ~isempty(ensemble.populations(end).models(jx).fixedExch)
            fixedExchs = ensemble.populations(end).models(jx).fixedExch;
        else
            fixedExchs = [];
        end
        
        mcaResults.enzNames = rxnNames;
        model = ensemble.populations(end).models(particleIdx(jx));
        if ix == 1
            xopt = ones(freeVars,1);
            xconst = ones(numel(ensemble.metsFixed), 1);
            vref = feval(kineticFxn,xopt,xconst,model,fixedExchs,Sred,kinInactRxns,subunits,0);
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
        xconst  = ones(numel(ensemble.metsFixed), numel(ix_mets));
        
        % Simulate flux for metabolite perturbation
        simFlux = feval(kineticFxn,xstep,xconst,model,fixedExchs,Sred,kinInactRxns,subunits,0);
        
        % Compute elasticiy matrices
        E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))'; % equivalent to imag(simFlux)./1.0e-10 ? 
        
        % Normalize elasticity matrices
        E_x_nor   = diag(vref.^-1)*E_x_abs*diag(xref);        
              
        % Compute control coefficients
        C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
        C_x       = diag(xref.^-1)*C_x_abs*diag(vref);
        C_v       = eye(numel(vref)) + E_x_nor*C_x;
        C_v(vref==0,:) = 0;                             % Make zero reactions with zero flux
        
        % Save control coefficients only if the result is accurate
        if all(abs(sum(C_x,2))<1e-5)
            if saveResMatrices
                mcaResults.xControl{ix}    = [mcaResults.xControl{ix}; C_x];
            end
            mcaResults.xControlAvg{ix} = mcaResults.xControlAvg{ix} + C_x;
            mcaResults.xcounter{ix}    = mcaResults.xcounter{ix} + 1;
        end
        if all(abs(sum(C_v(vref~=0,:),2))-1<1e-5)
            if saveResMatrices
                mcaResults.vControl{ix}    = [mcaResults.vControl{ix}; C_v];
                mcaResults.E_x_nor{ix}     = [mcaResults.E_x_nor{ix}; E_x_nor];
            end
            
            mcaResults.vControlAvg{ix} = mcaResults.vControlAvg{ix} + C_v;
            mcaResults.vcounter{ix}    = mcaResults.vcounter{ix} + 1;
        end
    end
    
    % Determine expectancy for control coefficients
    mcaResults.xControlAvg{ix} = mcaResults.xControlAvg{ix}/mcaResults.xcounter{ix};
    mcaResults.vControlAvg{ix} = mcaResults.vControlAvg{ix}/mcaResults.vcounter{ix};
    
end

end

