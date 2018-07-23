function [ensemble, models] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx)
%--------------------------------------------------------------------------
% Function used to calculate the reversibilities for each reaction
%
%-------------------------Marta Matos 2018---------------------------------

% Initialize reverTemp
ensemble.reverTemp = cell(size(ensemble.populations(1).probParams(strucIdx).rxnParams));
for activRxnIdx = 1:numel(ensemble.populations(1).probParams(strucIdx).rxnParams)
    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = zeros(size(ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx}'));
end

for activRxnIdx = 1:numel(ensemble.kinActRxns)     
		
    % Case 1: Diffusion and Exchanges
    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
            
        promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
        revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        alphaReversibility = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;
        
        numSamples = 1;
        thinning = 1;
        tol = 1e-10;
        
        % If the reaction is promiscuous
        if size(promiscRxnsList) > 0 
            
            % If the promiscuous reaction is the first on the list, calculate reversibilites
            if ensemble.kinActRxns(activRxnIdx) == promiscRxnsList(1)                                               
                gibbsTemp = cell2mat(ensemble.gibbsTemp(promiscRxnsList))';
                reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
                
                reverTemp2 = ones(1, size(revMatrix, 2));
                reverTemp2(sum(revMatrix) ~= 0) = reverTemp;
                
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
                
                % Back calculate normalized reversibilities
                randomRev = zeros(size(revMatrix));
                for rowI = 1:size(revMatrix,1)
                    randomRev(rowI, revMatrix(rowI,:) ~= 0) = log(reverTemp2(revMatrix(rowI,:)~=0));
                    randomRev(rowI, :) = randomRev(rowI, :) / (gibbsTemp(rowI)/RT);
                    assert(sum(randomRev(rowI, :)) < (1 + tol) && sum(randomRev(rowI, :)) > (1 - tol), ['Reversibilities do not sum up to 1, ', num2str(sum(randomRev(rowI, :)))]);
                end
                
                models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
            
            % Otherwise just copy the calculated reversibilities
            else                                                                                                     
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promiscRxnsList(1)};
                models(1).rxnParams(activRxnIdx).reversibilities = models(1).rxnParams(promiscRxnsList(1)).reversibilities;
            end
        
        % If the enzyme mechanism has multiple tracks, e.g. random
        elseif size(revMatrix,1)>1
            gibbsTemp = repelem(ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)}, size(revMatrix,1));
            reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
                
            reverTemp2 = ones(1, size(revMatrix, 2));
            reverTemp2(sum(revMatrix) ~= 0) = reverTemp;

            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
            
            % Back calculate normalized reversibilities
            randomRev = zeros(size(revMatrix));
            for rowI = 1:size(revMatrix,1)
                randomRev(rowI, revMatrix(rowI,:) ~= 0) = log(reverTemp2(revMatrix(rowI,:)~=0));
                randomRev(rowI, :) = randomRev(rowI, :) / (gibbsTemp(rowI)/RT);
                assert(sum(randomRev(rowI, :)) < (1 + tol) && sum(randomRev(rowI, :)) > (1 - tol), ['Reversibilities do not sum up to 1, ', num2str(sum(randomRev(rowI, :)))]);
            end

            models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
        
        % If the mechanism has only one track, e.g. ordered
        else 
            gibbsTemp = ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
            reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
            
            reverTemp2 = ones(1, size(revMatrix, 2));
            reverTemp2(revMatrix ~= 0) = reverTemp;

            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
                        
            % Back calculate normalized reversibilities
            randomRev = zeros(size(revMatrix));
            randomRev(revMatrix ~= 0) = log(reverTemp); 
            randomRev = randomRev / (gibbsTemp/RT);
            
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
            
            assert(sum(randomRev) < (1 + tol) && sum(randomRev) > (1 - tol), ['Reversibilities do not sum up to 1, ', num2str(sum(randomRev))]);
        end
      
    end
end

end

