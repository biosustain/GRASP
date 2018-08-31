function prevParams = getPopulationParamsMGK(ensemble,popIdx)
%--------------------------------------------------------------------------
% Assemble parameters from the previous population
%
% Inputs:       ensemble (structure), popIdx (double)
%
% Outputs:      preParams (structure array)
%--------------------- Pedro Saa 2017 -------------------------------------
% Initialize structure array
for strucIdx = 1:ensemble.numStruct
	prevParams(strucIdx).paramArray = [];											% This are the parameters transformed using ILRT
	prevParams(strucIdx).indexes    = [];											% Save parameter coordinates (indexes)
end

% Loop through all the structures in the model
for strucIdx = 1:ensemble.numStruct

	% Check that the structure is valid
	if (ensemble.populations(popIdx-1).structureWeights(strucIdx)==0)
		continue;

	% Proceed for this model structure
	else	
		aliveParticles = find(ensemble.populations(popIdx-1).strucIdx==strucIdx);
		prevIdx        = 1;
		
		% Loop through all the alive particles
		for jx = 1:numel(aliveParticles)
            
			% Initialize paramArray vector
			paramArray = [];
			
			% 1. Compute weights for pool factors (if any)
			if ~isempty(ensemble.poolConst)
                for ix = 1:numel(ensemble.poolConst)
                    poolFactorPrev      = ensemble.populations(popIdx-1).models(aliveParticles(jx)).poolFactor{ix};
                    poolFactorPrevTrans = map2uncons(poolFactorPrev);
                    paramArray          = [paramArray;poolFactorPrevTrans(:)];
                    if (jx == 1)
                        currIdx = prevIdx + numel(poolFactorPrevTrans) - 1;
                        prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
                        prevIdx = currIdx + 1;
                    end
                end
			end
	
			% 2. Compute weights for the gibbs factors
			if ~isempty(ensemble.thermoActive)
				gibbsFactorPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).gibbsFactor;
				gibbsFactorPrevTrans = log(gibbsFactorPrev./(1-gibbsFactorPrev));
				paramArray      = [paramArray;gibbsFactorPrevTrans(:)];
				if (jx == 1)
					currIdx = prevIdx + numel(gibbsFactorPrevTrans) - 1;
					prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
					prevIdx = currIdx + 1;
				end
			end
     
			% Loop thorugh all the reactions and initialize hyperparameters
			for activRxnIdx = 1:numel(ensemble.kinActRxns)
            
				% Skip non-enzymatic rxns
				if strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'diffusion')||...
					strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'freeExchange')||...
					strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'fixedExchange')||...
					strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'massAction')
					continue;
	
				% Compute importance weight for the rxn-specific parameters
                else					
				
					% 1. Compute weights of the reversibilities
					revMatrix    = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    revPrev      = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).reversibilities;
					revPrevTrans = map2uncons(revPrev);
					paramArray   = [paramArray;revPrevTrans(:)];
					if (jx == 1)
						currIdx = prevIdx + numel(revPrevTrans) - 1;
						prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
						prevIdx = currIdx + 1;
                    end
                    
                    % 2. Compute weights of the enzyme abundances
                    enzRPrev      = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).enzymeAbundances;
                    enzRPrevTrans = map2uncons(enzRPrev);
                    paramArray    = [paramArray;enzRPrevTrans(:)];
                    if (jx == 1)
                        currIdx = prevIdx + numel(enzRPrevTrans) - 1;
                        prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
                        prevIdx = currIdx + 1;
                    end
										
					% 3. Compute weights for branch factors (if any)
					Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
					if (size(Nelem,2)>1)						
						for ix = 1:size(Nelem,2)
							branchPrev(ix,1) = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).branchFactor(ix);                            
                            if (jx == 1)
                                currIdx = prevIdx + numel(branchPrev(ix)) - 1;
                                prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
                                prevIdx = currIdx + 1;
                            end
						end
						branchPrevTrans = log(branchPrev./(1-branchPrev));																						% Apply logit transformation
						paramArray      = [paramArray;branchPrevTrans(:)];
					end
																				
					% 4. Compute weights for mod elem flux (if any)
                    if (size(revMatrix,1)==1)&&any(revMatrix==0)
                        for ix = 1:sum(revMatrix==0)
                            modElemPrev(ix,1) = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
                            if (jx == 1)
                                currIdx = prevIdx + numel(modElemPrev) - 1;
                                prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
                                prevIdx = currIdx + 1;
                            end
                        end
                        modElemPrevTrans = log(modElemPrev./(1-modElemPrev));																					% Apply logit transformation
                        paramArray       = [paramArray;modElemPrevTrans(:)];                        
                    end

					% Check whether the current reaction is allosteric
					if ensemble.allosteric{strucIdx}(activRxnIdx)            
						allostericPrev      = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).allostericFactors;				% These parameters must be transformed
						allostericPrevTrans = log(allostericPrev./(1-allostericPrev));
						paramArray          = [paramArray;allostericPrevTrans(:)];
						if (jx == 1)
							currIdx = prevIdx + numel(allostericPrevTrans) - 1;
							prevParams(strucIdx).indexes = [prevParams(strucIdx).indexes;[prevIdx,currIdx]];
							prevIdx = currIdx + 1;
						end
					end
				end
			end
			
			% Save parameter array
			prevParams(strucIdx).paramArray = [prevParams(strucIdx).paramArray,paramArray];
		end			
	end       
end