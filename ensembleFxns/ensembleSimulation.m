function [fluxSS,fvalIter,exitFlag] = ensembleSimulation(model,ensemble,options)
% Testing ensemble of models for  perturbations - Pedro Saa Higuera 2015 UQ
% Description:  This function estimates the fluxes of an ensemble of models
%               for different perturbations
% Inputs:       ensemble of models, kineticFxn target, screenedModels,
%               reactions to be perturbed and perturbation levels
% Outputs:      relative differences of fluxes and steady-state fluxes

% Defining initial parameters
y0                = ones(sum(ensemble.metsSimulated),1);                         % Metabolite initial conditions
fluxKinetics      = str2func(ensemble.kineticFxn);                               % Specific kinetic function
[fluxNum,condNum] = size(ensemble.expressionLevels);
fluxSS            = zeros(fluxNum,condNum);
fvalIter          = zeros(1,condNum);
exitFlag          = zeros(1,condNum);
% Performing the perturbations to each screened model
for j = 1:condNum
    % For each model we solve for the perturbated reactions
    alphaTemp = ensemble.expressionLevels(:,j);    
    % Calculation of the normalized perturbed state
    [metSS,fval,exitflag] = fmincon(fluxKinetics,y0,[],[],[],[],1e-4*y0,1e3*y0,[],options,alphaTemp,model,ensemble,1);    
    fluxSS(:,j) = fluxKinetics(metSS,alphaTemp,model,ensemble,0);    
    fvalIter(j) = fval;
    exitFlag(j) = exitflag;
end