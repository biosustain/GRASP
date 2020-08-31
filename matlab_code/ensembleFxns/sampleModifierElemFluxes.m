function [models] = sampleModifierElemFluxes(ensemble, models, strucIdx)
% Function used to sample modifier elementary fluxes for each reaction.
%
% A modifier is any non-allosteric inhibitor or activator, and a modifier
% elementary flux is the flux of an elementary reaction
% in which a modifier is involved. 
% For instance, the flux of the elementary reaction in which a 
% non-allosteric inhibitor binds to the enzyme is a modifier elementary
% flux.
%
%
% USAGE:
%
%    models = sampleModifierElemFluxes(ensemble, models, strucIdx)
%
% INPUT:
%    ensemble (struct):	  model ensemble, see buildEnsemble for fields description
%    models (struct):     model, see initialSampler for fields description
%    strucIdx (int):      number of the model structure considered
%
% OUTPUT:
%    models (struct):	model structure with added reversibilities, see initialSampler for fields description
%
% .. Authors:
%       - Pedro Saa         2016 original code
%       - Marta Matos       2018 generalized it for promiscuous reactions
%                           and reactions with random mechanisms

for activRxnIdx = 1:numel(ensemble.kinActRxns)     
		
    % Case 1: Diffusion and Exchanges
    if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'fixedExchange')&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')

        promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};            
        revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        reverTemp = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};

        modifierElemFlux = [];
        
        % If the reaction is promiscuous and not the first on in the promiscRxnsList
        if size(promiscRxnsList,1) > 0 && ensemble.kinActRxns(activRxnIdx) ~= promiscRxnsList(1)
            modifierElemFlux = models(1).rxnParams(promiscRxnsList(1)).modiferElemFlux;
        
        % If the reaction is promiscuous and has inhibition steps
        elseif size(promiscRxnsList,1) > 0 && any(sum(revMatrix)==0)
            
            nElemFluxes = size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux, 1) ;
            modifierElemFlux = zeros(nElemFluxes, 2);
            
            for ix = 1:nElemFluxes
                aModifier              = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:));                    
                modifierElemFlux(ix,1) = aModifier(1)/sum(aModifier);       
                modifierElemFlux(ix,2) = modifierElemFlux(ix,1);
            end
        
        % If the reaction is not promiscuous and has inhibition steps
        % Marta: changed this to make it work for random mechanisms
        elseif ((size(revMatrix,1)==1) && any(revMatrix==0))  || ((size(revMatrix,1)>1) && any(sum(revMatrix)==0))
            modifierElemFlux = zeros(sum(reverTemp==1),2);
            for ix = 1:sum(reverTemp==1)
                aModifier              = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:));                    
                modifierElemFlux(ix,1) = aModifier(1)/sum(aModifier);
                modifierElemFlux(ix,2) = modifierElemFlux(ix,1);
            end
        end
        
        models(1).rxnParams(activRxnIdx).modiferElemFlux = modifierElemFlux;                       % save transpose of mod elem flux
    end
end

end

