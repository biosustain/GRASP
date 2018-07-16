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
        
        % if reaction is promiscuous and has no common intermediates and no
        % inhibitions
        if  size(alphaReversibility, 1) > 1  
            %disp('a');
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
        % if the reaction is promiscuous, and has an inhibition step
        %elseif size(promisc_rxns_list,1) > 0 && any(sum(revMatrix)==0)
        %    disp('d');
        %    disp(ensemble.rxns(ensemble.kinActRxns(activRxnIdx)));
        %    randomRev(sum(revMatrix)~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
        %    randomRev(sum(revMatrix)==0) = 0;
        elseif (size(revMatrix,1)==1)    
            %disp('b');
            models(1).rxnParams(activRxnIdx).reversibilities = randomRev';                                           % Save transpose
            randomRev(revMatrix~=0) = randomRev;                                                                     % Fill rev's for deadends with zeros
            randomRev(revMatrix==0) = 0;
        elseif (size(revMatrix,1)>1)                                                                              % Determine whether there are branches for the reversibility calculation
            %disp('c');
            lbRev     = randomRev;
            randomRev = computeBranchedReversibilities(revMatrix,lbRev);
            models(1).rxnParams(activRxnIdx).reversibilities = lbRev';                                               % Save transpose (note this has to be the lower bound!)
        end
        
        rTol = 1e-8;
        
        %assert(sum([sum(sumRandomRev < 1-tol) sum(sumRandomRev > 1+tol)])==0, 'Reversibilities do not sum up to one');

        % If the reaction is promiscuous and do not share enzyme intermediates,
        % except for free enzyme
        if size(promisc_rxns_list,1) > 0 && size(alphaReversibility, 1) > 1

            for rxnI = 1:size(promisc_rxns_list,2)                
                if promisc_rxns_list(rxnI) == ensemble.kinActRxns(activRxnIdx)
                    reverTemp = (revMatrix(rxnI,:)').*exp(randomRev(:,rxnI)*gibbsTemp/RT); 
                    %reverTemp = exp(randomRev(:,rxnI)*gibbsTemp/RT); 
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:, rxnI) = reverTemp;
                    
                    logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
                    logRev(logRev==-Inf) = 0;
                    %disp(rTol * gibbsTemp/RT);
                    assert(abs(sum(logRev(:,rxnI).*revMatrix(rxnI,:)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
                else
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(:,rxnI) = ensemble.reverTemp{promisc_rxns_list(rxnI)}(:,rxnI);
                end
            end

            if ensemble.kinActRxns(activRxnIdx) == promisc_rxns_list(numel(promisc_rxns_list))
                for rxnI = 1:(size(promisc_rxns_list,2)-1)
                    ensemble.reverTemp{promisc_rxns_list(rxnI)} = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};
                end
            end
            
        % If the reaction is promiscuous and shares enzyme intermediates    
        elseif size(promisc_rxns_list,1) > 0 && size(alphaReversibility, 1) == 1
            if ensemble.kinActRxns(activRxnIdx) == promisc_rxns_list(1)
                
                for rxnI = 1:size(revMatrix,1) 
                    %disp('in');
                    gibbsTemp = ensemble.gibbsTemp{promisc_rxns_list(rxnI)};
                    gibbsRT = exp(randomRev*gibbsTemp/RT);
                    ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)}(revMatrix(rxnI,:)'==1) = gibbsRT(revMatrix(rxnI,:)'==1); %exp(randomRev*gibbsTemp/RT);
                    
                    logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
                    logRev(logRev==-Inf) = 0;
                    %disp(sum(logRev.*revMatrix(rxnI,:)'));
                    %disp(gibbsTemp/RT);
                    %disp((abs(sum(logRev.*revMatrix(rxnI,:)') - gibbsTemp/RT)));
                    %disp(rTol * gibbsTemp/RT);
                    assert(abs(sum(logRev.*revMatrix(rxnI,:)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
                    
                end
                %logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
                %disp(logRev);
                %disp(sum(logRev.*revMatrix(1,:)'));
                %disp(gibbsTemp/RT);
                %assert(abs(sum(logRev.*revMatrix(1,:)') - gibbsTemp/RT) < 1e-10, 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
            else
                %disp('bla');
                ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = ensemble.reverTemp{promisc_rxns_list(1)};
            end
                
 
        else
            ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)} = exp(randomRev*gibbsTemp/RT);  % Convert to the proper units for later calculation  
            
            logRev = log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)});
            %disp(logRev);
            %disp(revMatrix(1, :));
            
            for rowI = 1:size(revMatrix, 1)
                %disp(rTol * gibbsTemp/RT);
                assert(abs(sum(logRev.*revMatrix(rowI, :)') - gibbsTemp/RT) < abs(rTol * gibbsTemp/RT), 'Sum of log-reversibilities does not add up to the Gibbs energy'); 
            end
        end
        %disp('---');
        %disp(gibbsTemp/RT);
        %for colI = 1:size(revMatrix(colI,:)
        %disp(sum(log(ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)})*revMatrix));
        %disp('***');
    end
end

end

