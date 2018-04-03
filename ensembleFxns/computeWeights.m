function ensemble = computeWeights(ensemble,popIdx)
%--------------------------------------------------------------------------
% Compute (unnormalized) weights based on current kernel and population
%
% Inputs:       ensemble (structure) ,population indicator (double)
%
% Outputs:      ensemble (structure)
%--------------------- Pedro Saa 2016 -------------------------------------
% Initialize weights to compute
weights = zeros(sum(ensemble.populations(popIdx).weights==0),1);    % Initialize weights
replenishedParticles = [];
aliveParticles       = [];

% Loop through all the strucutres
for strucIdx = 1:ensemble.numStruct

    % Find the particle indexes that belong to jx
    inactiveParticles = find((ensemble.populations(popIdx).strucIdx==strucIdx).*(ensemble.populations(popIdx).weights==0));
        
    % Skip calculations if there is no particles to replenish belonging to the current model structure
    if isempty(inactiveParticles); continue; end;
    
    % Determine alive particles for weight computation
    activeParticles = find((ensemble.populations(popIdx).strucIdx==strucIdx).*(ensemble.populations(popIdx).weights~=0));
    
    % Save replenished and alive mparticles
    replenishedParticles = [replenishedParticles;inactiveParticles];
    aliveParticles       = [aliveParticles;activeParticles];
    
    % Loop through particles belonging to jx
    for kx = 1:numel(inactiveParticles)
        
        % Initialize prior probability (model probability first)
        priorProb = (1/ensemble.numStruct);
        
        %% Compute prior probability of the current proposed particle
        % 1. Compute weights for pool factors (if any)
        if ~isempty(ensemble.poolConst)            
            poolTemp   = ensemble.populations(popIdx).models(inactiveParticles(kx)).poolFactor;
            alphaPrior = ensemble.populations(1).probParams(strucIdx).alphaPoolFactor;
            priorProb  = priorProb*probDirichlet(alphaPrior,poolTemp,1);
        end                

        % Loop thorugh all the reactions and initialize hyperparameters
        for activRxnIdx = 1:numel(ensemble.kinActRxns)

            % Initialize always except in the case of difusion mechanisms
            if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'difusion')&&~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')
                
                % 2. Compute weights for the Gibbs factors (only if necessary)
                if (ensemble.gibbsRanges(activRxnIdx,1)~=ensemble.gibbsRanges(activRxnIdx,2))
                    gibbsFactorTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).gibbsFactor;
                    betaPrior       = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaGibbsFactor;
                    priorProb       = priorProb*probDirichlet(betaPrior,[gibbsFactorTemp,1-gibbsFactorTemp],1);
                end
                if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction'); continue; end;
				
                % 3. Compute weights of the enzyme abundances
                enzRTemp   = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).enzymeAbundances;
                alphaPrior = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;
                priorProb  = priorProb*probDirichlet(alphaPrior,enzRTemp,1);

                % 4. Compute weights of the reversibilities
                revMatrix    = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
                revPrev      = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).reversibilities;
                revPrevValid = revPrev(revPrev~=0);                                                                     % remove fixed reversibilities
                alphaPrior   = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;
                priorProb    = priorProb*probDirichlet(alphaPrior,revPrevValid,1);

                % 5. Compute weights for branch factors (if any)
                Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
                if (size(Nelem,2)>1)
                    for ix = 1:size(Nelem,2)
                        branchTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).branchFactor(ix);
                        betaPrior  = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:);                        
                        priorProb  = priorProb*probDirichlet(betaPrior,[branchTemp,1-branchTemp],1);
                    end
                end

                % 6. Compute weights for mod elem flux (if any)
                if (size(revMatrix,1)==1)&&any(revMatrix==0)
                    for ix = 1:sum(revMatrix==0)
                        modElemTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
                        betaPrior   = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:);
                        priorProb   = priorProb*probDirichlet(betaPrior,[modElemTemp,1-modElemTemp],1);
                    end
                end

                % Check whether the current reaction is allosteric
                if ensemble.allosteric{strucIdx}(activRxnIdx)
                    
                    % 7. Compute weights for enzyme abundances (T state)
                    enzTemp    = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).enzymeAbundancesT;
                    alphaPrior = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundancesT;
                    priorProb  = priorProb*probDirichlet(alphaPrior,enzTemp,1);

                    % 8. Compute weights for relative activity
                    relActTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).relActivity;
                    betaPrior  = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaRelActivity;
                    priorProb  = priorProb*probDirichlet(betaPrior,[relActTemp,1-relActTemp],1);

                    % 9. Compute weights for the effector abundances (if any)
                    if ~isempty(ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)}) || ~isempty(ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)})
                        effBoundFractTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).effBoundFractions;
                        for ix = 1:max(size(effBoundFractTemp))
                            betaPrior = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaEffBoundFractions(ix,:);
                            priorProb = priorProb*probDirichlet(betaPrior,[effBoundFractTemp(ix),1-effBoundFractTemp(ix)],1);

                        end
                    end
                end
            end
        end

        % Initialize denominator weight of the current inactive particle
        denomWeight = 0;
        
        %% Compute transition probability for the current proposed particle
        for jx = 1:numel(activeParticles)
            
            % Initialize weight for the particular alive particle
            tempWeight = 1;
            
            % 1. Compute weights for pool factors (if any)
            if ~isempty(ensemble.poolConst)
                poolFactorPrev = ensemble.populations(popIdx).models(activeParticles(jx)).poolFactor;
                KpoolFactor    = ensemble.populations(popIdx).probParams(strucIdx).KpoolFactor;
                poolTemp       = ensemble.populations(popIdx).models(inactiveParticles(kx)).poolFactor;                
                tempWeight     = tempWeight*probDirichlet(KpoolFactor*poolFactorPrev+1,poolTemp,1);
            end
     
            % Loop thorugh all the reactions and initialize hyperparameters
            for activRxnIdx = 1:numel(ensemble.kinActRxns)
            
                % Initialize always except in the case of difusion mechanisms
                if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'difusion')&&~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')
                    
                    % 2. Compute weights for the Gibbs factors
                    if (ensemble.gibbsRanges(activRxnIdx,1)~=ensemble.gibbsRanges(activRxnIdx,2))
                        gibbsFactorPrev = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).gibbsFactor;
                        KgibbsFactor    = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KgibbsFactor;
                        gibbsTemp       = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).gibbsFactor;
                        tempWeight      = tempWeight*probDirichlet(KgibbsFactor*[gibbsFactorPrev,1-gibbsFactorPrev]+1,[gibbsTemp,1-gibbsTemp],1);
                    end
					if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction'); continue; end;
					
                    % 3. Compute weights of the enzyme abundances
                    enzRPrev   = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).enzymeAbundances;
                    KenzymesR  = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KenzymeAbundances;
                    enzRTemp   = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).enzymeAbundances;
                    tempWeight = tempWeight*probDirichlet(KenzymesR*enzRPrev+1,enzRTemp,1);

                    % 4. Compute weights of the reversibilities
                    revMatrix        = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    revPrev          = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).reversibilities;
                    Kreversibilities = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).Kreversibilities;
                    revTemp          = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).reversibilities;
                    tempWeight       = tempWeight*probDirichlet(Kreversibilities*revPrev+1,revTemp,1);

                    % 5. Compute weights for branch factors (if any)
                    Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    if (size(Nelem,2)>1)
                        for ix = 1:size(Nelem,2)
                            branchPrev    = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).branchFactor(ix);
                            KbranchFactor = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KbranchFactor;
                            branchTemp    = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).branchFactor(ix);
                            tempWeight    = tempWeight*probDirichlet(KbranchFactor*[branchPrev,1-branchPrev]+1,[branchTemp,1-branchTemp],1);
                        end
                    end

                    % 6. Compute weights for mod elem flux (if any)
                    if (size(revMatrix,1)==1)&&any(revMatrix==0)
                        for ix = 1:sum(revMatrix==0)
                            modElemPrev  = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
                            KmodElemTemp = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KmodiferElemFlux;
                            elemTemp     = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
                            tempWeight   = tempWeight*probDirichlet(KmodElemTemp*[modElemPrev,1-modElemPrev]+1,[elemTemp,1-elemTemp],1);
                        end
                    end

                    % Check whether the current reaction is allosteric
                    if ensemble.allosteric{strucIdx}(activRxnIdx)
                        
                        % 7. Compute weights for enzyme abundances (T state)
                        enzTPrev   = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).enzymeAbundancesT;
                        KenzymesT  = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KenzymeAbundancesT;
                        enzTemp    = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).enzymeAbundancesT;
                        tempWeight = tempWeight*probDirichlet(KenzymesT*enzTPrev+1,enzTemp,1);

                        % 8. Compute weights for relative activity
                        relActPrev   = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).relActivity;
                        KrelActivity = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KrelActivity;
                        relTemp      = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).relActivity;
                        tempWeight   = tempWeight*probDirichlet(KrelActivity*[relActPrev,1-relActPrev]+1,[relTemp,1-relTemp],1);

                        % 9. Compute weights for the effector abundances (if any)
                        if ~isempty(ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)}) || ~isempty(ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)})
                            effBoundFractTemp = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).effBoundFractions;
                            
                            for ix = 1:max(size(effBoundFractTemp))
                                effBoundFractPrev  = ensemble.populations(popIdx).models(activeParticles(jx)).rxnParams(activRxnIdx).effBoundFractions(ix);
                                KeffBoundFractions = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KeffBoundFractions;
                                effTemp            = ensemble.populations(popIdx).models(inactiveParticles(kx)).rxnParams(activRxnIdx).effBoundFractions(ix);
                                tempWeight         = tempWeight*probDirichlet(KeffBoundFractions*[effBoundFractPrev,1-effBoundFractPrev]+1,[effTemp,1-effTemp],1);                    
                            end
                        end
                    end                    
                end
            end

            % Add contribution of the current active to the inactive particle
            denomWeight = denomWeight + (ensemble.populations(popIdx).weights(activeParticles(jx))*tempWeight);            
        end

        % Compute model probability for the proposed particle
        modelProb = (ensemble.populations(popIdx-1).structureWeights(strucIdx));
        
        % Compute weight of current inactive particle
        weights(inactiveParticles(kx)) = priorProb/(denomWeight*modelProb);
    end
end

% Assign new weights and perform final normalization
disp(['Weights calculated for replenished particles: ',num2str(numel(replenishedParticles))]);
ensemble.populations(popIdx).weights(replenishedParticles) = (numel(replenishedParticles)/ensemble.numParticles)*weights/sum(weights);                                      % Normalize replenished weights
ensemble.populations(popIdx).weights(aliveParticles)       = (numel(aliveParticles)/ensemble.numParticles)*ensemble.populations(popIdx).weights(aliveParticles);            % Normalize weights of active particles