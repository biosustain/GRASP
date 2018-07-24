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
        tol = 1e-6;
        rTol = 1e-10;
        
        % If the reaction is promiscuous
        if size(promiscRxnsList) > 0 
         
            % If the promiscuous reaction is not the first in the list,just copy the reversibilities from the first
            if ensemble.kinActRxns(activRxnIdx) ~= promiscRxnsList(1)                                                                                                      
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promiscRxnsList(1)};
                models(1).rxnParams(activRxnIdx).reversibilities = models(1).rxnParams(promiscRxnsList(1)).reversibilities;
            
            % Otherwise calculate reversibilities
            else
                
                % If promiscuous reactions have steps in common
                if any(sum(revMatrix) >1)

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

                % If the promiscuous reactions do not have any steps in common
                else    
                    % Calculate normalized reversibilities
                    revTemp            = randg(alphaReversibility');
                    randomRev = zeros(size(revTemp));
                    for i = 1:size(revTemp, 2)
                        randomRev(:, i) = revTemp(:,i) / sum(revTemp(:,i));
                    end
                    randomRev = sum(randomRev, 2);

                    randomRev(sum(revMatrix)~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
                    randomRev(sum(revMatrix)==0) = 0;
                    models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose

                    % Calculate reversibilities
                    gibbsTemp = cell2mat(ensemble.gibbsTemp(promiscRxnsList))';
                    
                    reverTemp = zeros(size(revMatrix'));
                    for rxnI = 1:size(revMatrix, 1)
                        reverTemp(:, rxnI) = exp(randomRev*gibbsTemp(rxnI)/RT).*revMatrix(rxnI,:)';  % Convert to the proper units for later calculation   

                        logRev = log(reverTemp);
                        logRev(logRev==-Inf) = 0;
                        assert(abs(sum(logRev(:,rxnI).*revMatrix(rxnI,:)') - gibbsTemp(rxnI)/RT) < abs(rTol * gibbsTemp(rxnI)/RT), ['Sum of log-reversibilities does not add up to the Gibbs energy', num2str(abs(sum(logRev(:,rxnI).*revMatrix(rxnI,:)') - gibbsTemp(rxnI)/RT))]); 
                    end
                    
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = sum(reverTemp,2);
                end
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
        elseif size(revMatrix,1) == 1 
            % Sample normalized reversibilities
            revTemp = randg(alphaReversibility');
            randomRev = zeros(size(revTemp));
            randomRev = revTemp/sum(revTemp);
                                    
            randomRev(revMatrix~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
            randomRev(revMatrix==0) = 0;
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
            
            % Calculate reversibilities
            gibbsTemp = ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = exp(randomRev*gibbsTemp/RT);  % Convert to the proper units for later calculation  
            
            logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
            logRev(logRev==-Inf) = 0;
            assert(abs(sum(logRev.*revMatrix') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), ['Sum of log-reversibilities does not add up to the Gibbs energy', num2str(abs(sum(logRev.*revMatrix') - gibbsTemp/RT))]); 

        else
            disp('Oooops, something went wrong');
            exit
        end
      
    end
end

end

