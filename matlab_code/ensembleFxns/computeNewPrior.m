function computeNewPrior(ensemble)
%--------------------------------------------------------------------------
% Compute MLE of prob params for next prior
%
% Inputs:       ensemble (structure) , popIdx (double) , idxValidParticles
%               (cell array)
%
% Outputs:      ensemble with MLE values computed 9structure)
%--------------------- Pedro Saa 2016 -------------------------------------
% Initialize optimization parameters for FMINCON
optionsMLE = optimset('Algorithm','sqp','Display','off','TolFun',1e-12,'TolX',1e-10,'MaxFunEvals',2e4,'MaxIter',2e4);
c = '%';

% Loop through all the structures and rxns in the model
for strucIdx = 1:ensemble.numStruct
    disp(['Adjusting prior based on current posterior: Model ',num2str(strucIdx),' in progress...']);      
    
    % Extract alive particles from the current model structure
    idxParticleStructure = find((ensemble.populations(end).strucIdx==strucIdx));
    
    % Skip structure with no valid particles
    if isempty(idxParticleStructure); continue; end;
    
    % I. Extract weights from previous population    
    weightsIdxStruc = ensemble.populations(end).weights(idxParticleStructure);
    
    % Write new script
    fid = fopen(['Adjusted_prior',num2str(strucIdx),'.txt'],'w'); 
    
    % 1. Compute MLE for the pool factors (if any)
    if ~isempty(ensemble.poolConst)
        
        % Retrieve information from the previous population
        poolFactorPrev = zeros(length(idxParticleStructure),size(ensemble.populations(end).probParams(strucIdx).alphaPoolFactor,2));
        for ix = 1:length(idxParticleStructure)
            poolFactorPrev(ix,:) = ensemble.populations(end).models(idxParticleStructure(ix)).poolFactor';
        end
        
        % Solve for the Dir parameters
        alpha0 = ensemble.populations(end).probParams(strucIdx).alphaPoolFactor;
        try
            alphaPool = fmincon(@probDirichlet,alpha0,[],[],[],[],zeros(size(alpha0)),[],[],optionsMLE,poolFactorPrev,2,weightsIdxStruc');             % Solve the constrained problem
        catch                                                                             % If the above doesn't work, solve the unconstrained problem using a derivative-free algorithm
            alphaPool = fminsearch(@probDirichlet,log(alpha0),optionsMLE,poolFactorPrev,3,weightsIdxStruc');                                           % Solve the unconstrained problem
            alphaPool = exp(alphaPool);
        end
        fprintf(fid,'alphaPool %d\n',alphaPool);
    end
    
    % II. Loop thorugh all the reactions and initialize hyperparameters
    for activRxnIdx = 1:length(ensemble.kinActRxns)
        
        % Initialize always except in the case of difusion mechanisms
        if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'difusion')&&~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')            
            fprintf(fid,'%s\n',ensemble.rxns{ensemble.kinActRxns(activRxnIdx)});
            
            % 0. Compute MLE for gibbs factors and pool factors (if any)
            gibbsFactorPrev = zeros(length(idxParticleStructure),1);
            for ix = 1:length(idxParticleStructure)
                gibbsFactorPrev(ix,1) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).gibbsFactor;
            end           
            
            % Solve Beta parameters for the Gibbs factors
            beta0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaGibbsFactor;
            try
                betaGibbs = fmincon(@probDirichlet,beta0,[],[],[],[],zeros(size(beta0)),[],[],optionsMLE,[gibbsFactorPrev,1-gibbsFactorPrev],2,weightsIdxStruc');       % Solve the constrained problem
            catch                                                                                          % If the above doesn't work, solve the unconstrained problem using a derivative-free algorithm
                betaGibbs = fminsearch(@probDirichlet,log(beta0),optionsMLE,[gibbsFactorPrev,1-gibbsFactorPrev],3,weightsIdxStruc');                                    % Solve the unconstrained problem
                betaGibbs = exp(betaGibbs);
            end            
            fprintf(fid,'betaGibbs %d\n',betaGibbs);
            
            % If the mechanism is mass action, continue to the next iteration
            if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction'); continue; end;
            
            % 1. Retrieve EnzR information from the previous population
            enzRPrev = zeros(length(idxParticleStructure),size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances,2));
            for ix = 1:length(idxParticleStructure)
                enzRPrev(ix,:) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).enzymeAbundances;
            end
            
            % Solve Dir parameters for the enzyme abundances
            alpha0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;%
            try
                alphaEnzR = fmincon(@probDirichlet,alpha0,[],[],[],[],zeros(size(alpha0)),[],[],optionsMLE,enzRPrev,2,weightsIdxStruc');                                % Solve the constrained problem
            catch                                                                                          % If the above doesn't work, solve the unconstrained problem using a derivative-free algorithm
                alphaEnzR = fminsearch(@probDirichlet,log(alpha0),optionsMLE,enzRPrev,3,weightsIdxStruc');                                                              % Solve the unconstrained problem
                alphaEnzR = exp(alphaEnzR);
            end
            fprintf(fid,'alphaEnzR %d\n',alphaEnzR);            
            
            % 2. Retrieve Rev information from the previous population
            revPrev = zeros(length(idxParticleStructure),size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities,2));
            for ix = 1:length(idxParticleStructure)
                revTemp       = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).reversibilities;
                revPrev(ix,:) = revTemp(revTemp~=0);
            end
            
            % Solve Dir parameters for the reversibilities
            alpha0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;
            try
                alphaRev = fmincon(@probDirichlet,alpha0,[],[],[],[],zeros(size(alpha0)),[],[],optionsMLE,revPrev,2,weightsIdxStruc');                                  % Solve the constrained problem
            catch                                                                                         % If the above doesn't work, solve the unconstrained problem using a derivative-free algorithm
                alphaRev = fminsearch(@probDirichlet,log(alpha0),optionsMLE,revPrev,3,weightsIdxStruc');                                                                % Solve the unconstrained problem
                alphaRev = exp(alphaRev);
            end
            fprintf(fid,'alphaRev %d\n',alphaRev);           
            
            % 3. Retrieve branch information from the previous population (if any)
            Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
            
            % If there is branch information available
            if (size(Nelem,2)>1)
                branchPrev = zeros(length(idxParticleStructure),size(Nelem,2));
                for ix = 1:length(idxParticleStructure)
                    branchPrev(ix,:) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).branchFactor;
                end
                
                % Solve Beta parameters for the branch factors
                for ix = 1:size(branchPrev,2)
                    beta0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:);
                    try
                        betaBranch = fmincon(@probDirichlet,beta0,[],[],[],[],zeros(size(beta0)),[],[],optionsMLE,[branchPrev(:,ix),1-branchPrev(:,ix)],2,weightsIdxStruc');
                    catch
                        betaBranch = fminsearch(@probDirichlet,log(beta0),optionsMLE,[branchPrev(:,ix),1-branchPrev(:,ix)],3,weightsIdxStruc');
                        betaBranch = exp(betaBranch);
                    end
                    fprintf(fid,'betaBranch %d\n',betaBranch);
                end
            end
            
            % 4. Retrieve information about modifiers flux (if any)            
            if (size(ensemble.revMatrix{activRxnIdx,strucIdx},1)==1)&&any(ensemble.revMatrix{activRxnIdx,strucIdx}==0)
                modElemFluxPrev = zeros(length(idxParticleStructure),sum(revTemp==0));
                for ix = 1:length(idxParticleStructure)
                    modElemFluxPrev(ix,:) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).modiferElemFlux;
                end
                
                % Solve Beta parameters for the branch factors
                for ix = 1:sum(revTemp==0)
                    beta0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:);
                    try
                        betaModifier = fmincon(@probDirichlet,beta0,[],[],[],[],zeros(size(beta0)),[],[],optionsMLE,[modElemFluxPrev(:,ix),1-modElemFluxPrev(:,ix)],2,weightsIdxStruc');
                    catch
                        betaModifier = fminsearch(@probDirichlet,log(beta0),optionsMLE,[modElemFluxPrev(:,ix),1-modElemFluxPrev(:,ix)],3,weightsIdxStruc');
                        betaModifier = exp(betaModifier);
                    end
                    fprintf(fid,'betaModifier %d\n',betaModifier);
                end
            end
            
            % Check whether the current reaction is allosteric
            if ensemble.allosteric{strucIdx}(activRxnIdx)
                
                % 5/6. Retrieve EnzT and rel. activity information from the previous population
                enzTPrev   = zeros(length(idxParticleStructure),size(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundancesT,2));
                relPrevAct = zeros(length(idxParticleStructure),1);
                for ix = 1:length(idxParticleStructure)
                    enzTPrev(ix,:) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).enzymeAbundancesT;
                    relPrevAct(ix) = ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).relActivity;
                end
                
                % 5. Solve Dir parameters for the enzyme abundances
                alpha0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundancesT;
                try
                    alphaEnzT = fmincon(@probDirichlet,alpha0,[],[],[],[],zeros(size(alpha0)),[],[],optionsMLE,enzTPrev,2,weightsIdxStruc');
                catch
                    alphaEnzT = fminsearch(@probDirichlet,log(alpha0),optionsMLE,enzTPrev,3,weightsIdxStruc');
                    alphaEnzT = exp(alphaEnzT);
                end
                fprintf(fid,'alphaEnzT %d\n',alphaEnzT);
                
                % 6. Solve for the relative activity ~ Beta(a,b)
                beta0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaRelActivity;
                try
                    betaRelActivity = fmincon(@probDirichlet,beta0,[],[],[],[],zeros(size(beta0)),[],[],optionsMLE,[relPrevAct,1-relPrevAct],2,weightsIdxStruc');
                catch
                    betaRelActivity = fminsearch(@probDirichlet,log(beta0),optionsMLE,[relPrevAct,1-relPrevAct],3,weightsIdxStruc');
                    betaRelActivity = exp(betaRelActivity);
                end
                fprintf(fid,'betaRelActivity %d\n',betaRelActivity);
                
                % 7. Retrieve information about allosteric effectors (if any)
                if ~isempty(ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)}) ||...
                        ~isempty(ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)})
                    
                    % Extract information from the previous population
                    effBoundFractPrev = zeros(length(idxParticleStructure),size(ensemble.populations(end).models(idxParticleStructure(1)).rxnParams(activRxnIdx).effBoundFractions,2));
                    for ix = 1:length(idxParticleStructure)
                        effBoundFractPrev(ix,:) =  ensemble.populations(end).models(idxParticleStructure(ix)).rxnParams(activRxnIdx).effBoundFractions;
                    end
                    
                    % Computation beta parameters
                    for ix = 1:size(effBoundFractPrev,2)
                        beta0 = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaEffBoundFractions(ix,:);
                        try
                            betaEffector = fmincon(@probDirichlet,beta0,[],[],[],[],zeros(size(beta0)),[],[],optionsMLE,[effBoundFractPrev(:,ix),1-effBoundFractPrev(:,ix)],2,weightsIdxStruc');
                        catch
                            betaEffector = fminsearch(@probDirichlet,log(beta0),optionsMLE,[effBoundFractPrev(:,ix),1-effBoundFractPrev(:,ix)],3,weightsIdxStruc');
                            betaEffector = exp(betaEffector);
                        end
                        fprintf(fid,'betaEffector %d\n',betaEffector);
                    end
                end
            end
        end
    end
    fclose(fid);    
end
disp('MLE computations performed.');