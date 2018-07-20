function [ensemble, models] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx)


% Initialize reverTemp
ensemble.reverTemp = cell(size(ensemble.populations(1).probParams(strucIdx).rxnParams));
for activRxnIdx = 1:numel(ensemble.populations(1).probParams(strucIdx).rxnParams)
    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = zeros(size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities'));
end

for activRxnIdx = 1:numel(ensemble.kinActRxns)     
    %disp(ensemble.rxns(ensemble.kinActRxns(activRxnIdx)));
		
    % Case 1: Diffusion and Exchanges
    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
            
        
        
        promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
        revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        alphaReversibility = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;
        
        numSamples = 1;
        thinning = 1;
        tol = 1e-12;
        %DGr1    = -10;                  % [kJ/mol]
        
        % if the reaction is promiscuous
        if size(promiscRxnsList) > 0 
        
            if ensemble.kinActRxns(activRxnIdx) == promiscRxnsList(1)                                               % if the promiscuous reaction is the first on the list, calculate reversibilites
                gibbsTemp = cell2mat(ensemble.gibbsTemp(promiscRxnsList))';
                reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
                
                reverTemp2 = ones(1, size(revMatrix, 2));
                reverTemp2(sum(revMatrix) ~= 0) = reverTemp;
                
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
                
                % back calculate reversibilities
                randomRev = zeros(size(revMatrix));
                for rowI = 1:size(revMatrix,1)
                    randomRev(rowI, revMatrix(rowI,:) ~= 0) = log(reverTemp2(revMatrix(rowI,:)~=0));
                    randomRev(rowI, :) = randomRev(rowI, :) / (gibbsTemp(rowI)/RT);
                    assert(sum(randomRev(rowI, :)) < (1 + tol) && sum(randomRev(rowI, :)) > (1 - tol), 'Reversibilities do not sum up to 1');
                end
                
                models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
  
            else                                                                                                     % otherwise just copy the calculated reversibilities
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promiscRxnsList(1)};
                models(1).rxnParams(activRxnIdx).reversibilities = models(1).rxnParams(promiscRxnsList(1)).reversibilities;
            end
        
        % if mechanism with multiple tracks, e.g. random
        elseif size(revMatrix,1)>1
            gibbsTemp = repelem(ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)}, size(revMatrix,1));
            reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
                
            reverTemp2 = ones(1, size(revMatrix, 2));
            reverTemp2(sum(revMatrix) ~= 0) = reverTemp;

            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
            
            % back calculate reversibilities
            randomRev = zeros(size(revMatrix));
            for rowI = 1:size(revMatrix,1)
                randomRev(rowI, revMatrix(rowI,:) ~= 0) = log(reverTemp2(revMatrix(rowI,:)~=0));
                randomRev(rowI, :) = randomRev(rowI, :) / (gibbsTemp(rowI)/RT);
                assert(sum(randomRev(rowI, :)) < (1 + tol) && sum(randomRev(rowI, :)) > (1 - tol), 'Reversibilities do not sum up to 1');
            end

            models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
        
        else 
            gibbsTemp = ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
            reverTemp = generalRevSampling(alphaReversibility, gibbsTemp/RT, numSamples, thinning);
            
            reverTemp2 = ones(1, size(revMatrix, 2));
            reverTemp2(revMatrix ~= 0) = reverTemp;

            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = reverTemp2';
                        
            % back calculate reversibilities
            randomRev = zeros(size(revMatrix));
            randomRev(revMatrix ~= 0) = log(reverTemp); 
            randomRev = randomRev / (gibbsTemp/RT);
            
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev;
            
            assert(sum(randomRev) < (1 + tol) && sum(randomRev) > (1 - tol), 'Reversibilities do not sum up to 1');
        end
        
        
        

        %

        %revTemp            = randg(alphaReversibility');

        %randomRev = zeros(size(revTemp));
        %if size(alphaReversibility, 1) > 1 && size(promiscRxnsList, 1) > 0
        %    for i = 1:size(revTemp, 2)
        %        randomRev(:, i) = revTemp(:,i) / sum(revTemp(:,i));
        %    end
        %else
        %    randomRev          = revTemp/sum(revTemp);
        %end
        
        % if reaction is promiscuous and has no common intermediates and no
        % inhibitions
%         if  size(alphaReversibility, 1) > 1  
%             models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
%         elseif (size(revMatrix,1)==1)    
%             models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
%             randomRev(revMatrix~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
%             randomRev(revMatrix==0) = 0;
%         elseif (size(revMatrix,1)>1)        
%             lbRev     = randomRev;
%             randomRev = computeBranchedReversibilities(revMatrix,lbRev);
%             models(1).rxnParams(activRxnIdx).reversibilities = lbRev';                                               % Save transpose (note this has to be the lower bound!)
%         end
%         
%         rTol = 1e-8;
%         
%         %assert(sum([sum(sumRandomRev < 1-tol) sum(sumRandomRev > 1+tol)])==0, 'Reversibilities do not sum up to one');
% 
%         % If the reaction is promiscuous and do not share enzyme intermediates,
%         % except for free enzyme
%         if size(promiscRxnsList,1) > 0 && size(alphaReversibility, 1) > 1
% 
%             for rxnI = 1:size(promiscRxnsList,2)                
%                 
%                 if promiscRxnsList(rxnI) == ensemble.kinActRxns(activRxnIdx)
%                     
%                     reverTemp = (revMatrix(rxnI,:)').*exp(randomRev(:,rxnI)*gibbsTemp/RT); 
%                     
%                     ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:, rxnI) = reverTemp;
%                     
%                     logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
%                     logRev(logRev==-Inf) = 0;
%                     assert(abs(sum(logRev(:,rxnI).*revMatrix(rxnI,:)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
%                
%                 else
%                     ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:,rxnI) = ensemble.reverTemp{promiscRxnsList(rxnI)}(:,rxnI);
%                 end
%             end
% 
%             
%             
%         % If the reaction is promiscuous and shares enzyme intermediates    
%         elseif size(promiscRxnsList,1) > 0 && size(alphaReversibility, 1) == 1
%             if ensemble.kinActRxns(activRxnIdx) == promiscRxnsList(1)
%                 
%                 for rxnI = 1:size(revMatrix,1) 
%                     gibbsTemp = ensemble.gibbsTemp{promiscRxnsList(rxnI)};
%                     gibbsRT = exp(randomRev*gibbsTemp/RT);
%                     ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(revMatrix(rxnI,:)'==1) = gibbsRT(revMatrix(rxnI,:)'==1);
%                     
%                     logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
%                     logRev(logRev==-Inf) = 0;
% 
%                     assert(abs(sum(logRev.*revMatrix(rxnI,:)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
%                     
%                 end
%                 
%             else
%                 ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promiscRxnsList(1)};
%             end
%                 
%  
%         else
%             ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = exp(randomRev*gibbsTemp/RT);  % Convert to the proper units for later calculation        
%             logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
% 
%             for rowI = 1:size(revMatrix, 1)
%                 assert(abs(sum(logRev.*revMatrix(rowI, :)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
%             end
%         end

    end
end

end

