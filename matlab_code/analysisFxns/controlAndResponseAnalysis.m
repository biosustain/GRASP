function mcaResults = controlAndResponseAnalysis(ensemble,saveResMatrices,conditionI,strucIdx)
% Calculates both control coefficients and response coefficients.
% Response coefficients here answer the question: what happens if one 
% increases a given enzyme's concentration?
%
% This is mostly important when modeling promiscuous enzymes, where 
% increasing the enzyme's concentration doesn't necessarily lead to an 
% equal increase in the flux of the reactions it catalyzes.
% Response and control coefficients are the same when enzymes are
% independent and an increase in enzyme concentration leads to a
% proportional increase in the reaction flux.
%
% Response coefficients are calculated as 
%
% .. math::
%           C_E^J = C_v^J\Pi
%
% where :math:`C_v^J` is the flux control coefficient matrix and 
% :math:`\Pi` is the parameter elasticity matrix.
%
% The implementation is based on
%
%  - https://doi.org/10.1002/9780470475935.ch3, section on MCA
%  - https://doi.org/10.1111/j.1432-1033.1990.tb15329.x
%  - https://doi.org/10.1111/j.1432-1033.1990.tb15330.x
%
% For all models in the given ensemble it returns the average flux and 
% concentration control coefficients across the ensemble as well as the
% average enzyme and concentration response coefficients.
%
% It can also return the control and response coefficient matrices for
% each model in the ensemble if *saveResMatrices* is set to true. But keep 
% in mind that for large models this can use a lot of RAM and CPU.
%
%
% USAGE:
%
%    mcaResults = controlAndResponseAnalysis(ensemble, saveResMatrices, strucIdx)
%
% INPUT:
%    ensemble (struct):           model ensemble, see buildEnsemble for fields description
%    saveResMatrices (logical):   whether or not to save the elasticity and control coefficient matrices for all models
%
% OPTIONAL INPUT:
%    conditionI (int):  condition for which the user wants for run MCA when having multiple coniditions, by
%                       default MCA is run for all conditions
%    strucIdx (int):    number of the model structure considered
%
% OUTPUT:
%    mcaResults (struct):	control and response analysis results
%
%               * xControlAvg (*cell*)   : average concentration controls coefficient for the model ensemble
%               * vControlAvg (*cell*)   : average flux control coefficients for the model ensemble
%               * eResponseAvg (*cell*)	 : average enzyme response coefficient for the model ensemble
%               * xResponseAvg (*cell*)	 : average concentration response coefficients for the model ensemble
%               * xcounter (*cell*)      : number of models in the average concentration control coefficient calculation
%               * vcounter (*cell*)      : number of models in the average flux control coefficient calculation
%               * xRcounter (*cell*)     : number of models in the average concentration response coefficient calculation
%               * eRcounter (*cell*)     : number of models in the average enzyme response coefficient calculation
%               * xControl (*cell*)      : concentration control coefficient matrix for each model
%               * vControl (*cell*)      : flux control coefficient matrix for each model
%               * xResponse (*cell*)     : concentration response coefficient matrix for each model
%               * eResponse (*cell*)     : enzyme response coefficient matrix for each model
%
% .. Authors:
%       - Pedro Saa     2018 original code
%       - Marta Matos   2019 refactored code and added response analysis


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
numModels   = numel(particleIdx);

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
        mcaResults.xControl{ix}     = [];
        mcaResults.vControl{ix}     = [];
        mcaResults.eResponse{ix}    = [];
        mcaResults.xResponse{ix}    = [];
    end
    
    mcaResults.xControlAvg{ix}  = 0;
    mcaResults.vControlAvg{ix}  = 0;
    mcaResults.eResponseAvg{ix} = 0;
    mcaResults.xResponseAvg{ix} = 0;
    mcaResults.xcounter{ix}     = 0;
    mcaResults.vcounter{ix}     = 0;
    mcaResults.eRcounter{ix}    = 0;
    mcaResults.xRcounter{ix}    = 0;
    
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
        
        % Define step length to perturb enzyme concentrations
        hstep_e      = hstep*Eref;
        xmets        = repmat(xref,1,numel(Eref));
        xenz         = repmat(Eref,1,numel(Eref));
        xenzImag     = 1i*diag(hstep_e);
        origXenzImag = xenzImag;
        
        for rxnI=1:size(ensemble.promiscuity{ix},2)
            if size(ensemble.promiscuity{ix}{rxnI},2) ~= 0
                sumVector = zeros(numel(Eref),1);
                for promiscuousRxnI=1:size(ensemble.promiscuity{ix}{rxnI},2)
                     sumVector = sumVector + origXenzImag(:,ensemble.promiscuity{ix}{rxnI}(promiscuousRxnI));
                end
                xenzImag(:,rxnI) = sumVector;
            end
        end
        
        xenz    = xenz + xenzImag;
        estep   = [xmets;xenz];
        xconst  = ones(numel(ensemble.metsFixed), size(estep,2));
        
        % Simulate flux for enzyme perturbation
        simFluxEnz = feval(kineticFxn,estep,xconst,model,fixedExchs,Sred,kinInactRxns,subunits,0);
        
        % Compute elasticiy matrices
        E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))'; % equivalent to imag(simFlux)./1.0e-10 ? 
        E_pi_abs = -(imag(simFluxEnz')./hstep_e(:,ones(1,numFluxes)))';% equivalent to imag(simFluxEnz)./1.0e-10 ? 
        
        % Normalize elasticity matrices
        E_x_nor   = diag(vref.^-1)*E_x_abs*diag(xref);
        E_pi_nor  = diag(vref.^-1)*E_pi_abs*diag(Eref);
        
        % Delete repeated columns
        columnsToDelete = [];
        for rxnI=1:size(ensemble.promiscuity{ix},2)
            if size(ensemble.promiscuity{ix}{rxnI},2) ~= 0
                for promiscuousRxnI=2:size(ensemble.promiscuity{ix}{rxnI},2)
                    columnsToDelete = [columnsToDelete, ensemble.promiscuity{ix}{rxnI}(promiscuousRxnI)];
                end
            end
        end
        columnsToDelete = sort(unique(columnsToDelete), 'descend');
                   
        for colI=columnsToDelete
            E_pi_nor(:, colI)           = [];
            mcaResults.enzNames(colI,:) = [];
        end
                
        % Compute control coefficients
        C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
        C_x       = diag(xref.^-1)*C_x_abs*diag(vref);
        C_v       = eye(numel(vref)) + E_x_nor*C_x;
        C_v(vref==0,:) = 0;  
        R_e       = C_v * E_pi_nor;
        R_x       = C_x * E_pi_nor;
        
        % Save control coefficients only if the result is accurate
        if all(abs(sum(C_x,2))<1e-5)
            if saveResMatrices
                mcaResults.xControl{ix}     = [mcaResults.xControl{ix}; C_x];
                mcaResults.xResponse{ix}    = [mcaResults.xResponse{ix}; R_x];
            end
            
            mcaResults.xControlAvg{ix}  = mcaResults.xControlAvg{ix} + C_x;
            mcaResults.xcounter{ix}     = mcaResults.xcounter{ix} + 1;
            mcaResults.xResponseAvg{ix} = mcaResults.xResponseAvg{ix} + R_x;
            mcaResults.xRcounter{ix}    = mcaResults.xRcounter{ix} + 1;
        end
        if all(abs(sum(C_v(vref~=0,:),2))-1<1e-5)
            if saveResMatrices
                mcaResults.vControl{ix}     = [mcaResults.vControl{ix}; C_v];
                mcaResults.eResponse{ix}    = [mcaResults.eResponse{ix}; R_e];
            end
            
            mcaResults.vControlAvg{ix}  = mcaResults.vControlAvg{ix} + C_v;
            mcaResults.vcounter{ix}     = mcaResults.vcounter{ix} + 1;
            mcaResults.eResponseAvg{ix} = mcaResults.eResponseAvg{ix} + R_e;
            mcaResults.eRcounter{ix}    = mcaResults.eRcounter{ix} + 1;
        end
    end
    
    % Determine expectancy for control coefficients
    mcaResults.xControlAvg{ix} = mcaResults.xControlAvg{ix}/mcaResults.xcounter{ix};
    mcaResults.vControlAvg{ix} = mcaResults.vControlAvg{ix}/mcaResults.vcounter{ix};
    mcaResults.eResponseAvg{ix} = mcaResults.eResponseAvg{ix}/mcaResults.eRcounter{ix};
    mcaResults.xResponseAvg{ix} = mcaResults.xResponseAvg{ix}/mcaResults.xRcounter{ix};
end

end
