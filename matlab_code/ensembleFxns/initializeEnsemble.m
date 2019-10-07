function ensemble = initializeEnsemble(ensemble,popIdx,verbose)
% Initialize ensemble structure.
%
%
% USAGE:
%
%    ensemble = initializeEnsemble(ensemble, popIdx, verbose)
%
% INPUT:
%    ensemble (struct):  model ensemble
%    popIdx (int):       population ID
%
% OPTIONAL INPUT:
%    verbose (int):      verbosity level
%                                  
% OUTPUT:
%    ensemble (struct):   initialized model ensemble data structure
%
%               * description (*char*)            : model name basically
%               * sampler (*char*)                : specifies the sampling mode: ORACLE or rejection [TODO Pedro check if this is correct]
%               * solver (*char*)                 : which solver to use for the rejection sampler
%               * numConditions (*int*)           : how many experimental conditions    
%               * numStruct (*int*)               : how many model structures
%               * numParticles (*int*)            : how many models to be sampled   
%               * parallel (*int*)                : whether or not to parallelize sampling process
%               * numCores (*int*)                : if *parallel* is true, how many cores to use
%               * alphaAlive (*double*)           : [TODO Pedro]
%               * tolerance (*double vector*)     : [TODO Pedro]
%               * S (*int matrix*)                : stoichiometric matrix as defined in the input file
%               * rxns (*char cell*)              : reaction IDs
%               * rxnNames (*char cell*)          : reaction names
%               * exchRxns (*int vector*)         : exchange reactions, marked as transport in rxns sheet
%               * activeRxns (*int vector*)       : list with reactions marked as modeled
%               * isoenzymes (*cell*)             : list with isoenzymes
%               * uniqueIso (*cell*)              : [TODO Nick]
%               * mets (*char cell*)              : metabolite IDs
%               * metNames (*char cell*)          : metabolite names
%               * rxnMets (*char cell*)           : [TODO Pedro]
%               * metsBalanced (*int vector*)     : [TODO Pedro]  
%               * metsSimulated (*int vector*)    : [TODO Pedro]  
%               * metsFixed (*int vector*)        : which metabolites concentrations are defined as fixed (constant)
%               * measuredMets (*int vector*)     : which metabolite concentrations were actually measured (vs. estimated)
%               * Sred (*int matrix*)             : reduced stoichiometric matrix, includes only balanced metabolites and active reactions
%               * measRates (*double matrix*)     : measured reaction fluxes means
%               * measRatesStd (*double matrix*)  : measured reaction fluxes standard deviations
%               * splitRatios (*vector*)          : [TODO Pedro]  
%               * poolConst (*vector*)            : [TODO Pedro]  
%               * ineqThermoConst (*vector*)      : [TODO Pedro]      
%               * expFluxes (*double vector*)     : [TODO Pedro]  
%               * expFluxesStd (*double vector*)  : [TODO Pedro]     
%               * fluxRef (*double vector*)       : reference reaction fluxes means
%               * fluxRefStd (*double vector*)    : reference reaction fluxes standard deviations
%               * freeFluxes (*int vector*)       : [TODO Pedro]  
%               * simWeights (*double vector*)    : [TODO Pedro]  
%               * Sthermo (*int matrix*)          : stoichiometric matrix used for thermodynamics, excludes exchange reactions
%               * gibbsRanges (*double matrix*)   : thermodynamically feasible ranges for Gibbs energies  
%               * metRanges (*double matrix*)     : thermodynamically feasible ranges for metabolite concentrations
%               * G0Ranges (*double matrix*)      : thermodynamically feasible ranges for standard Gibbs energies
%               * metsDataMin (*double vector*)   : minimum value for scaled metabolite concentrations   
%               * metsDataMax (*double vector*)   : maximum value for scaled metabolite concentrations  
%               * metsDataMean (*double vector*)  : mean value for scaled metabolite concentrations   
%               * prevPrior (*cell*)              : [TODO Pedro] 
%               * prevPriorInfo (*cell*)          : [TODO Pedro]     
%               * allosteric (*cell*)             : which reactions are allosterically regulated
%               * subunits (*int cell*)           : number of enzyme subunits for each reaction
%               * rxnMechanisms (*char cell*)     : reaction mechanisms    
%               * extremePathways (*int cell*)    : extreme pathways for the given reaction mechanism
%               * inhibitors (*char cell*)        : reaction inhibitors
%               * activators (*char cell*)        : reaction activators 
%               * negEffectors (*char cell*)      : allosteric inhibitors
%               * posEffectors (*char cell*)      : allosteric activators   
%               * subOrder (*char cell*)          : substrate binding order
%               * prodOrder (*char cell*)         : product release order
%               * promiscuity (*int cell*)        : promiscuous reactions  
%               * kinActRxns (*int vector*)       : kinetically active reactions, includes all reactions with mechanism other than fixedExchange 
%               * prodDataMin (*double vector*)   : minimum value for scaled enzyme concentrations
%               * prodDataMax (*double vector*)   : maximum value for scaled enzyme concentrations  
%               * prodDataMean (*double vector*)  : mean value for scaled enzyme concentrations  
%               * kinInactRxns (*int vector*)     : kinetically inactive reactions, basically reactions with a fixedExchange mechanism 
%               * fixedExch (*double matrix*)     : fixed exchange reactions
%               * kineticFxn (*char cell*)        : name of kinetic function used to build the model with all rate laws
%               * metLists (*char cell*)          : list of metabolites (as defined in patterns) for each reaction
%               * revMatrix (*int matrix*)        : [TODO Pedro]
%               * forwardFlux (*int cell*)        : [TODO Pedro]  
%               * Nelem (*int cell*)              : [TODO Pedro]
%               * rxnMetLinks (*char cell*)       : [TODO Pedro]  
%               * freeVars (*char cell*)          : [TODO Pedro]
%               * metsActive (*int vector*)       : [TODO Pedro]
%               * eigThreshold (*double*)         : threshold for the real part of the jacobian eigenvalues, if there is any higher than the threshold the model is discarded
%               * thermoActive (*int vector*)     : [TODO Pedro]
%               * populations (*struct*)          : [TODO Pedro]
%                       * probParams (*struct*)
%                               * muGibbsFactor (*double vector*)     : [TODO Pedro]
%                               * sigmaGibbsFactor (*double vector*)  : [TODO Pedro]
%                               * rxnParams (*struct*)                : reaction parameters
%                                       * alphaEnzymeAbundances (*int vector*)  : [TODO Pedro]
%                                       * alphaReversibilities (*int vector*)   : [TODO Pedro]
%                                       * betaModifierElemeFlux (*int vector*)  : [TODO Pedro]
%                                       * betaBranchFactor (*int vector*)       : [TODO Pedro]
%                        * weights (*double vector*)  : [TODO Pedro]
%                        * models (*struct*)          : models in the ensemble
%                                * poolFactor (*double vector*)  : [TODO Pedro]
%                                * gibbsFactor (*double vector*) : [TODO Pedro]
%                                * rxnParams (*struct*)          : reaction parameters
%                                        * reversibilities (*double vector*)   : sampled elementary reaction reversibilities
%                                        * enzymeAbundances (*double vector*)  : sampled enzyme intermediates abundances
%                                        * branchFactor (*double vector*)      : [TODO Pedro]
%                                        * modifierElemeFlux (*double vector*) : [TODO Pedro]
%                                        * allostericFactors (*double vector*) : [TODO Pedro]
%                                        * L (*double vector*)                 : [TODO Pedro]
%                                        * KnegEff (*double vector*)           : [TODO Pedro]
%                                        * KposEff (*double vector*)           : [TODO Pedro]
%               * replenished particles (*int*)   : [TODO Pedro]
%
% .. Authors:
%       - Pedro Saa         2016 original code
%       - Marta Matos       2018 generalized it for promiscuous reactions 
%                           and random mechanisms
%       - Nicholas Cowie	2019 added extreme pathways


if (nargin<3); verbose = 0; end

% Determine the type of sampler selected (rejection is the only sampler currently supported)
if strcmpi(ensemble.sampler,'ORACLE')||strcmpi(ensemble.sampler,'rejection')||strcmpi(ensemble.sampler,'SMC_rejection')||strcmpi(ensemble.sampler,'SMC')||strcmpi(ensemble.sampler,'MCMC-SMC')

	% Determine thermo-active reactions
	ensemble.thermoActive = find((ensemble.gibbsRanges(:,1)~=-1e2)&(ensemble.gibbsRanges(:,1)~=ensemble.gibbsRanges(:,2)));

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
                    extremePathways = ensemble.extremePathways{strucIdx}{activRxnIdx};
                    if (size(extremePathways,1)>1)
                        ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor = ones(size(extremePathways,1),1);
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
        if any(ensemble.prevPrior); ensemble = extractPrior(ensemble); end
    end

    % Preallocate memory for sampling parameters
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

if verbose; disp('Ensemble successfully initialized.'); end
