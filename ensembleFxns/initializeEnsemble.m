function ensemble = initializeEnsemble(ensemble,popIdx,verbose)
%--------------------------------------------------------------------------
% Initialize ensemble structure
%
% Inputs:   (ensemble)    ensemble (structure)
%
% Outputs:  (ensemble)    initialized ensemble (structure)
%------------------------Pedro Saa 2016, Marta Matos 2018------------------
if (nargin<3); verbose = 0; end;

% Determine the type of sampler selected
if (strcmpi(ensemble.sampler,'rejection')||strcmpi(ensemble.sampler,'SMC'))
    
	% Determine thermo-active reactions
	ensemble.thermoActive = find((ensemble.gibbsRanges(:,1)~=-1e3)&(ensemble.gibbsRanges(:,1)~=ensemble.gibbsRanges(:,2)));
	
    % Find independent variables for the any of the sampling algorithm (only the first time)
    if (popIdx==1)       
	   	   
       % Initialize the remaining parameters based on the strucIdx and rxnIdx
        for strucIdx = 1:ensemble.numStruct           
            
            % Check whether there are any pool factors
            if ~isempty(ensemble.poolConst)
                for ix = 1:numel(ensemble.poolConst)
                    Aeq = ensemble.poolConst{ix}(1:numel(ensemble.metsActive));
                    ensemble.populations(1).probParams(strucIdx).alphaPoolFactor{ix} = ones(1,sum((Aeq~=0)));
                end
            end
			
			% Initialize hyper-parameters for the Gibbs free energy of reaction ~ MVN(0,sigma)
			ensemble.populations(1).probParams(strucIdx).muGibbsFactor    = zeros(numel(ensemble.thermoActive),1);
            ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor = eye(numel(ensemble.thermoActive));
			
            % Loop thorugh all the reactions and initialize hyperparameters
            for activRxnIdx = 1:numel(ensemble.kinActRxns)
                
                % Initialize always except in the case of exchange/difussion/mass action mechanisms
                if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')&&...
					~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')&&...
                    ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'fixedExchange')&&...
					~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
                                                            
                    % Initialize prob parameters of enzyme abundances ~ Dir(1)
                    forwardFlux = ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances = ones(1,max(forwardFlux(:)));
                    
                    % Initialize prob parameters of reversibilities (check that there is no dead-end reversibility) ~ Dir(1)
                    revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    
                    promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
                    % For promiscuous reactions set alphaReversibilities = revMatrix after removing columns with all zeroes 
                    if size(promiscRxnsList,1) > 0  
                        
                        % Remove zero columns from revMatrix
                        revMatrixTemp = revMatrix;
                        revMatrixTemp( :, ~any(revMatrixTemp,1) ) = [];
                        ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities = revMatrixTemp;
                        
                        % Initialize modifiers (if any) ~ numModifiers x Beta(1)
                        if any(sum(revMatrix)==0)
                            nRows = (sum(sum(revMatrix)==0));
                            ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux = ones(nRows,2);
                        end
                        
                    % If the reaction is not promiscuous
                    elseif (size(revMatrix,1)==1)
                        ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities = ones(1,sum(revMatrix));
                        
                        % Initialize modifiers (if any) ~ numModifiers x Beta(1)
                        if any(revMatrix==0)
                            ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux = ones(sum(revMatrix==0),2);
                        end
                        
                    % For branched mechanisms do as follows
                    % Marta: changing this to work with random mechanisms with inhibitions - hopefully it won't break
                    elseif (size(revMatrix,1)>1)
                        revMatrixTemp = revMatrix;
                        revMatrixTemp( :, ~any(revMatrixTemp,1) ) = [];
                        ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities = ones(1,size(revMatrixTemp,2));
                        
                        % Initialize modifiers (if any) ~ numModifiers x Beta(1)
                        if any(sum(revMatrix)==0)
                            nRows = (sum(sum(revMatrix)==0));
                            ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux = ones(nRows,2);
                        end
                        
                    end
                    
                    % Initialize branch factors ~ Beta(1)
                    Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
                    if (size(Nelem,2)>1)
                        ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor = ones(size(Nelem,2),2);                        
                    end
                    
                    % Initialize paramters for allosteric rxns
                    if ensemble.allosteric{strucIdx}(activRxnIdx)                                                
                        
                        % Initialize prob parameters for effectors (order: L, posEff & negEff)
                        posEffectors = ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
                        negEffectors = ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
						numEffectors = max(size(posEffectors)) + max(size(negEffectors));					
						ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).muAllostericParams    = zeros(numEffectors+1,1);
						ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams = eye(numEffectors+1);
                    end
                end
            end
        end
        
        % Set replenished particle number
        ensemble.replenishedParticles(popIdx) = ensemble.numParticles;       
        
        % Set unnormalized weights
        ensemble.populations(1).weights = ones(ensemble.replenishedParticles(popIdx),1); 
        
        % Check whether there are previous prior information
        if any(ensemble.prevPrior); ensemble = extractPrior(ensemble); end;
    end
    
    % Preallocate memory for sampled parameters
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).poolFactor(numel(ensemble.kinActRxns))                  = 0;    % pool constraint factor ~ Dir(1)
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).gibbsFactor(numel(ensemble.thermoActive))               = 0;    % gibbs free energy of reaction ~ MVN(0,1)
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).reversibilities   = 0;    % microscopic reversibilites, assumed equal for the active (R) and tense (T) states ~ Dir(1)
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).enzymeAbundances  = 0;    % enzyme abudances in the R state ~ Dir(1)    
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).branchFactor      = 0;    % branch factor ~ B(1) (size equal to size(Nelem,2)>1)
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).modiferElemFlux   = 0;    % modifier elem flux (for non-allosteric inhibitors and/or activators) ~ B(1)
	ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).allostericFactors = 0;    % allosteric factors used in the calculation of L, KposEff & KnegEff
        
    % Preallocate memory for computed parameters
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).L                 = 0;    % allosteric constant
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).KnegEff           = 0;    % allosteric binding constant activator
    ensemble.populations(popIdx).models(ensemble.replenishedParticles(popIdx)).rxnParams(numel(ensemble.kinActRxns)).KposEff           = 0;    % allosteric binding constant inhibitor
else
    disp('Sampler selected not found.'); return;
end
if verbose; disp('Ensemble successfully initialized.'); end;