function [models] = sampleModifierElemFluxes(ensemble, models, strucIdx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


for activRxnIdx = 1:numel(ensemble.kinActRxns)     
    %disp(ensemble.rxns(ensemble.kinActRxns(activRxnIdx)));
		
    % Case 1: Diffusion and Exchanges
    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')

        promisc_rxns_list = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};            
        gibbsTemp = ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
        revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        reverTemp = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};

        modifierElemFlux = [];

        if size(promisc_rxns_list,1) > 0 && ensemble.kinActRxns(activRxnIdx) ~= promisc_rxns_list(1)
            modifierElemFlux = models(1).rxnParams(promisc_rxns_list(1)).modiferElemFlux';

        elseif size(promisc_rxns_list,1) > 0 && any(sum(revMatrix)==0)
            
            nElemFluxes = size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:), 2) ;
            modifierElemFlux = zeros(nElemFluxes, 2);
            
            for ix = 1:nElemFluxes
                aModifier              = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:));                    
                modifierElemFlux(ix,1) = aModifier(1)/sum(aModifier);       
                modifierElemFlux(ix,2) = modifierElemFlux(ix,1);
            end
            modifierElemFlux = modifierElemFlux';

        else
            if ((size(revMatrix,1)==1) && any(revMatrix==0)) 
                modifierElemFlux = zeros(sum(reverTemp==1),1);
                for ix = 1:sum(reverTemp==1)
                    aModifier              = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:));                    
                    modifierElemFlux(ix,1) = aModifier(1)/sum(aModifier);                    
                end
            end
        end
        models(1).rxnParams(activRxnIdx).modiferElemFlux = modifierElemFlux';                       % save transpose of mod elem flux
    end
end

end

