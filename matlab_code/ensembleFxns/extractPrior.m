function ensemble = extractPrior(ensemble)
% Extract and assign new prior to the current ensemble .
%
%
% USAGE:
%
%    ensemble = extractPrior(ensemble)
%
% INPUT:
%    ensemble (struct):    model ensemble
%
% OUTPUT:
%    ensemble (struct):    model ensemble
%
% .. Authors:
%       - Pedro Saa     2016 original code

indexes = find(ensemble.prevPrior);
for strucIdx = indexes
    
    % Determine positions of the parameters of the reactions
    numParam = find(ensemble.prevPriorInfo{strucIdx,1}==0);
    
    % Find indexes of the reactions with priors
    for jx = 1:numel(ensemble.activeRxns)
        rxnPresent = strcmp(ensemble.rxns{ensemble.activeRxns(jx)},ensemble.prevPriorInfo{strucIdx,2});             % find reaction in the prior list
              
        % If there is a match: find position in the list ot extract parameters
        if any(rxnPresent)
            currentRxn = strcmp(ensemble.prevPriorInfo{strucIdx,2}(rxnPresent),ensemble.rxns(ensemble.kinActRxns));
            paramIdx   = find(strcmp(ensemble.prevPriorInfo{strucIdx,2},ensemble.rxns{ensemble.activeRxns(jx)}));
            endIndex   = numParam(numParam>paramIdx);
            
            % Figure out indexes for the parameters
            if ~isempty(endIndex)                
                tempParams = paramIdx+1:endIndex(1)-1;
            else
                tempParams = paramIdx+1:length(ensemble.prevPriorInfo{strucIdx,1});
            end
            
            % Extract parameter information and values
            paramIden = ensemble.prevPriorInfo{strucIdx,2}(tempParams);
            paramVals = ensemble.prevPriorInfo{strucIdx,1}(tempParams);
            
            % Start filling out prior parameter by parameter
            % 1. Gibbs free energy
            paramTemp = find(strcmp(paramIden,'betaGibbs'));
            if ~isempty(paramTemp)
                ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).betaGibbsFactor = paramVals(paramTemp)';
            end
            
            % 2. Enzyme abundances (active state)
            paramTemp = find(strcmp(paramIden,'alphaEnzR'));
            if ~isempty(paramTemp)
                ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).alphaEnzymeAbundances = paramVals(paramTemp)';
            end
            
            % 3. Reversibilities
            paramTemp = find(strcmp(paramIden,'alphaRev'));
            if ~isempty(paramTemp)
                ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).alphaReversibilities = paramVals(paramTemp)';
            end
            
            % 4. Branch factor
            paramTemp = find(strcmp(paramIden,'betaBranch'));            
            if ~isempty(paramTemp)
                lengthParam = length(paramTemp);
                for ix = 1:lengthParam/2
                    ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).betaBranchFactor(ix,:) = paramVals(paramTemp(2*ix-1:2*ix))';
                end
            end
            
            % 5. Modifiers
            paramTemp = find(strcmp(paramIden,'betaModifier'));
            if ~isempty(paramTemp)
                lengthParam = length(paramTemp);
                for ix = 1:lengthParam/2
                    ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).betaModiferElemFlux(ix,:) = paramVals(paramTemp(2*ix-1:2*ix))';
                end
            end
            
            % 6. Enzyme abundances (tense state)
            paramTemp = find(strcmp(paramIden,'alphaEnzT'));
            if ~isempty(paramTemp)
                ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).alphaEnzymeAbundancesT = paramVals(paramTemp)';
            end
                        
            % 7. Relative activity
            paramTemp = find(strcmp(paramIden,'betaRelActivity'));
            if ~isempty(paramTemp)
                ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).betaRelActivity = paramVals(paramTemp)';
            end
            
            % 8. Relative activity
            paramTemp = find(strcmp(paramIden,'betaEffector'));
            if ~isempty(paramTemp)
                lengthParam = length(paramTemp);
                for ix = 1:lengthParam/2
                    ensemble.populations(1).probParams(strucIdx).rxnParams(currentRxn).betaEffBoundFractions(ix,:) = paramVals(paramTemp(2*ix-1:2*ix))';
                end
            end            
        end       
    end
    
    % If there is information about the pool factors, incorporate that into the model prior
    poolIdxs = find(strcmp(ensemble.prevPriorInfo{strucIdx,2},'alphaPool'));
    if ~isempty(poolIdxs)
        ensemble.populations(1).probParams(strucIdx).poolParams = ensemble.prevPriorInfo{strucIdx,2}(poolIdxs(2:end));
    end
end
ensemble = rmfield(ensemble,'prevPriorInfo');
disp('Prior successfully loaded.');