function [ensemble] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx)


% Initialize reverTemp
ensemble.reverTemp = cell(size(ensemble.populations(1).probParams(strucIdx).rxnParams));
for activRxnIdx = 1:numel(ensemble.populations(1).probParams(strucIdx).rxnParams)
    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = zeros(size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities'));
end

for activRxnIdx = 1:numel(ensemble.kinActRxns)        
		
    % Case 1: Diffusion and Exchanges
    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
            
        gibbsTemp = ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
        
        promisc_rxns_list = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
        revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};

        alphaReversibility = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;

        revTemp            = randg(alphaReversibility');

        randomRev = zeros(size(revTemp));
        if size(alphaReversibility, 1) > 1 && size(promisc_rxns_list, 1) > 0
            for i = 1:size(revTemp, 2)
                randomRev(:, i) = revTemp(:,i) / sum(revTemp(:,i));
            end
        else
            randomRev          = revTemp/sum(revTemp);
        end

        if  size(alphaReversibility, 1) > 1  
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
        elseif (size(revMatrix,1)==1)    
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
            randomRev(revMatrix~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
            randomRev(revMatrix==0) = 0;
        elseif (size(revMatrix,1)>1)                                                                              % Determine whether there are branches for the reversibility calculation
            lbRev     = randomRev;
            randomRev = computeBranchedReversibilities(revMatrix,lbRev);
            models(1).rxnParams(activRxnIdx).reversibilities = lbRev';                                               % Save transpose (note this has to be the lower bound!)
        end

        % When promiscuous reactions do not share enzyme intermediates,
        % except for free enzyme
        if size(promisc_rxns_list,1) > 0 && size(alphaReversibility, 1) > 1

            for rxn_i = 1:size(promisc_rxns_list,2)                
                if promisc_rxns_list(rxn_i) == ensemble.kinActRxns(activRxnIdx)
                    reverTemp = (revMatrix(rxn_i,:)').*exp(randomRev(:,rxn_i)*gibbsTemp/RT); 
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:, rxn_i) = reverTemp;
                else
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:,rxn_i) = ensemble.reverTemp{promisc_rxns_list(rxn_i)}(:,rxn_i);
                end
            end

            if ensemble.kinActRxns(activRxnIdx) == promisc_rxns_list(numel(promisc_rxns_list))
                for rxn_i = 1:(size(promisc_rxns_list,2)-1)
                    ensemble.reverTemp{promisc_rxns_list(rxn_i)} = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};
                end
            end
        % When promiscuous reactions share enzyme intermediates    
        elseif size(promisc_rxns_list,1) > 0 && size(alphaReversibility, 1) == 1
            if ensemble.kinActRxns(activRxnIdx) == promisc_rxns_list(1)
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = exp(randomRev*gibbsTemp/RT);
            else
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promisc_rxns_list(1)};
            end
                
 
        else
            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = exp(randomRev*gibbsTemp/RT);  % Convert to the proper units for later calculation    
        end

    end
end

end

