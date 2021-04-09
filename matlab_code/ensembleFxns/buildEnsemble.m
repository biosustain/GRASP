function ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold)
% Samples a kinetic model ensemble which contains only valid models.
%
% Valid models are models where
%
%  - for all reactions the fluxes and respective Gibbs energies are 
%    compatible;
%  - the real part of the jacobian eigenvalues is lower than the defined 
%    threshold;
%  - the difference between the predicted flux and the reference flux
%    is negligible.
%
% It will sample models until there are *n* valid models, where n is
% specified in the input excel file as the number of particles.
% To avoid sampling forever, *maxNumberOfSamples* is defined and no more
% than *maxNumberOfSamples* models will be sampled, regardless of the number
% of valid models sampled.
% 
% Models with a real part of the jacobian eigenvalues greather than 
% *eigThreshold* are discarded.
%
% If the output directory doesn't exist it is created.
%
%
% USAGE:
%
%    ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples, eigThreshold)
%
% INPUT:
%    inputFile (char):            path to input file
%    outputFile (char):           path to output file
%    maxNumberOfSamples (int):    maximum number of models to be sampled
%    eigThreshold (double):       threshold for positive eigenvalues' real part
%
% OUTPUT:
%    ensemble (struct):    model ensemble
%
%               * description (*char*)            : model name basically
%               * sampler (*char*)                : specifies the sampling mode: GRASP or rejection
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
%               * freeVars (*char cell*)          : [TODO Pedro]
%               * metsActive (*int vector*)       : [TODO Pedro]
%               * eigThreshold (*double*)         : threshold for the real part of the jacobian eigenvalues, if there is any higher than the threshold the model is discarded
%               * thermoActive (*int vector*)     : [TODO Pedro]
%               * populations (*struct*)          : [TODO Pedro]
%
%                       * probParams (*struct*)      : [TODO Pedro]
%
%                               * muGibbsFactor (*double vector*)     : [TODO Pedro]
%                               * sigmaGibbsFactor (*double vector*)  : [TODO Pedro]
%                               * rxnParams (*struct*)                : reaction parameters
%
%                                       * alphaEnzymeAbundances (*int vector*)  : [TODO Pedro]
%                                       * alphaReversibilities (*int vector*)   : [TODO Pedro]
%                                       * betaModifierElemeFlux (*int vector*)  : [TODO Pedro]
%                                       * betaBranchFactor (*int vector*)       : [TODO Pedro]
%                       * weights (*double vector*)  : [TODO Pedro]
%                       * models (*struct*)          : models in the ensemble
%
%                                * poolFactor (*double vector*)  : [TODO Pedro]
%                                * gibbsFactor (*double vector*) : [TODO Pedro]
%                                * rxnParams (*struct*)          : reaction parameters
%
%                                        * reversibilities (*double vector*)   : sampled elementary reaction reversibilities
%                                        * enzymeAbundances (*double vector*)  : sampled enzyme intermediates abundances
%                                        * branchFactor (*double vector*)      : [TODO Pedro]
%                                        * modifierElemeFlux (*double vector*) : [TODO Pedro]
%                                        * kineticParams (*double vector*)     : reaction kinetic parameters
%                       * strucIdx (*double vector*)  : model structure ID
%                       * tolScore (*double vector*)  : [TODO Pedro]
%                       * xopt (*cell*)               : [TODO Pedro]
%                       * simFluxes (*cell*)          : [TODO Pedro]
%               * replenished particles (*int*)   : [TODO Pedro]
%
% .. Authors:
%       - Pedro Saa     2016 original code
%       - Marta Matos	2019 refactored code and added check for valid
%                       models

% 1. Load information
popIdx   = 1;
ensemble = loadEnsembleStructure(inputFile);
ensemble.eigThreshold = eigThreshold;


% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);


% 3. Sample fluxes and Gibbs energies
priorType = 'normal';
ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples,priorType);


% Create output folder if it doesn't exist yet
[filepath,~,~] = fileparts(outputFile);
if ~exist(filepath)
    mkdir(filepath);
end

disp('Running initial sampler.');

% Setup folder with temp files
try
    rmdir('temp','s');
    mkdir('temp');
catch
    mkdir('temp');
end

% Preallocate memory for the remaing fields in the ensemble structure
tolScore = zeros(ensemble.replenishedParticles(popIdx),ensemble.numConditions);
strucIdx = zeros(ensemble.replenishedParticles(popIdx),1);
xopt{ensemble.replenishedParticles(popIdx),1}      = [];
simFluxes{ensemble.replenishedParticles(popIdx),1} = [];


% Check whether the job is ran in parallel    
if ensemble.parallel

    sampleCount = 0;
    nValidModels = 0;

    while nValidModels < ensemble.numParticles

        if sampleCount == 0                                             % if no models have been sampled yet, sample ensemble.numParticles models
            nSamples = ceil(ensemble.numParticles * 1.3);
        else                                                            % else check how many more should be sampled based on the percentage of valid models
            nSamples = 1 / (nValidModels / sampleCount)*(sampleCount - nValidModels);
            nSamples = round(nSamples);

            if nSamples > maxNumberOfSamples - sampleCount
                nSamples = maxNumberOfSamples - sampleCount;
            end
        end

        parpool(ensemble.numCores);
        parfor ix = (sampleCount+1):(sampleCount+nSamples)   
            [validModelList(ix),models(ix),strucIdx(ix),xopt{ix},tolScore(ix,:),simFluxes{ix}] = initialSampler(ensemble, ix);
        end
        delete(gcp);

        sampleCount = sampleCount + nSamples;           
        nValidModels = size(models(validModelList), 2);

        if nValidModels == 0
            break
        end

    end

    if nValidModels > 0
        models = models(validModelList);
        strucIdx = strucIdx(validModelList);
        xopt = xopt{validModelList};
        tolScore = tolScore(validModelList);
        simFluxes = simFluxes{validModelList};            
    end

else
    sampleCount = 1;
    ix = 1;

    while ix <= ensemble.numParticles && sampleCount <= maxNumberOfSamples

        [isModelValid,model,strucIdxTemp,xoptTemp,tolScoreTemp,simFluxesTemp] = initialSampler(ensemble,sampleCount);

        if isModelValid
            models(ix) = model;
            strucIdx(ix) = strucIdxTemp;
            xopt{ix} = xoptTemp;
            tolScore(ix) = tolScoreTemp;
            simFluxes{ix} = simFluxesTemp;

            ix = ix + 1;
        end

        sampleCount = sampleCount + 1;
    end

    nValidModels = ix - 1;
    sampleCount = sampleCount -1;

end

disp([newline, newline, '*** Out of ',num2str(sampleCount), ' models, ', num2str(nValidModels), ' were valid ***', newline, newline]);


if nValidModels > 1
    % Append sampling results
    ensemble.populations(1).strucIdx  = strucIdx;                                                                           % model structures
    ensemble.populations(1).tolScore  = tolScore;                                                                           % tolerance score
    ensemble.populations(1).xopt      = xopt;                                                                               % optimal value found
    ensemble.populations(1).simFluxes = simFluxes;                                                                          % simulated fluxes
    ensemble.populations(1).models    = models;                                                                             % model particles
    ensemble = rmfield(ensemble, 'gibbsEnergies');
    ensemble = rmfield(ensemble, 'metConcRef');
    clearvars -except ensemble popIdx modelID outputFile
    save(outputFile);
else
    disp('No valid models were sampled.');
end

end

