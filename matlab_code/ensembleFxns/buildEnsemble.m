function ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold,testing)
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
%               * tolerance (*double*)            : max allowable flux discrepancy between simulations and flux data
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
%               * metsBalanced (*int vector*)     : Indices of mass-balanced metabolites
%               * metsFixed (*int vector*)        : which metabolites concentrations are defined as fixed (constant)
%               * Sred (*int matrix*)             : reduced stoichiometric matrix, includes only balanced metabolites and active reactions
%               * measRates (*double matrix*)     : measured reaction fluxes means
%               * measRatesStd (*double matrix*)  : measured reaction fluxes standard deviations
%               * poolConst (*vector*)            : pool conservation constraints 
%               * ineqThermoConst (*vector*)      : thermodynamic inequality constraints on log concentrations (ratios in linear space). Used in TMFA
%               * expFluxes (*double vector*)     : mean of experimental fluxes
%               * expFluxesStd (*double vector*)  : std of experimental fluxes    
%               * fluxRef (*double vector*)       : reference reaction fluxes means
%               * fluxRefStd (*double vector*)    : reference reaction fluxes standard deviations
%               * freeFluxes (*int vector*)       : indices of linearly independent fluxes   
%               * simWeights (*double vector*)    : discrepancy flux weights (by default experimental fluxes) 
%               * Sthermo (*int matrix*)          : stoichiometric matrix used for thermodynamics, excludes exchange reactions
%               * gibbsRanges (*double matrix*)   : thermodynamically feasible ranges for Gibbs energies  
%               * metRanges (*double matrix*)     : thermodynamically feasible ranges for metabolite concentrations
%               * G0Ranges (*double matrix*)      : thermodynamically feasible ranges for standard Gibbs energies
%               * metsDataMin (*double vector*)   : minimum value for scaled metabolite concentrations   
%               * metsDataMax (*double vector*)   : maximum value for scaled metabolite concentrations  
%               * metsDataMean (*double vector*)  : mean value for scaled metabolite concentrations   
%               * prevPrior (*cell*)              : previous prior information available for the reactions
%               * prevPriorInfo (*cell*)          : description of previous prior information  
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
%               * revMatrix (*int matrix*)        : reversibilty matrix for the kinetic mechanism
%               * forwardFlux (*int cell*)        : link matrix of forward elementary fluxes connecting enzyme intermediates
%               * Nelem (*double cell*)           : null space basis for the mass balanced elementary fluxes
%               * freeVars (*char cell*)          : free variables, includes both proteins and metabolites concentrations
%               * metsActive (*int vector*)       : indices of active (non-constant) metabolites
%               * eigThreshold (*double*)         : threshold for the real part of the jacobian eigenvalues, if there is any higher than the threshold the model is discarded
%               * thermoActive (*int vector*)     : thermodynamically active reactions indices (e.g. excludes fixed exchanges) 
%               * populations (*struct*)          : population of model ensembles
%
%                       * probParams (*struct*)      : probabilitis parameters of the ensemble
%
%                               * muGibbsFactor (*double vector*)     : [THIS SHOULD BE REMOVED FROM THE ENSEMBLE STRUCTURE AS WE USE THE HR SAMPLER]
%                               * sigmaGibbsFactor (*double vector*)  : [THIS SHOULD BE REMOVED FROM THE ENSEMBLE STRUCTURE AS WE USE THE HR SAMPLER]
%                               * rxnParams (*struct*)                : reaction parameters
%
%                                       * alphaEnzymeAbundances (*double vector*)  : prior hyper parameters of the Dirichlet distribution for sampling enzyme abundances
%                                       * alphaReversibilities (*double vector*)   : prior hyper parameters of the Dirichlet distribution for sampling reversibilities
%                                       * betaModifierElemeFlux (*double vector*)  : prior hyper parameters of the beta distribution for sampling elementary modifier flux
%                                       * betaBranchFactor (*double vector*)       : prior hyper parameters of the beta distribution for flux branching factors
%                       * weights (*double vector*)  : relative weight of the models
%                       * models (*struct*)          : models in the ensemble
%
%                                * poolFactor (*double vector*)  : proportions of conserved moieties at the reference state
%                                * rxnParams (*struct*)          : reaction parameters
%
%                                        * reversibilities (*double vector*)   : sampled elementary reaction reversibilities
%                                        * enzymeAbundances (*double vector*)  : sampled enzyme intermediates abundances
%                                        * branchFactor (*double vector*)      : random branching flux factor in random order mechanisms (one per branching point)
%                                        * modifierElemeFlux (*double vector*) : random flux through an elementary reaction involving a modifier (inhibitor or activator) 
%                                        * kineticParams (*double vector*)     : reaction kinetic parameters
%                       * strucIdx (*double vector*)  : model structure ID
%                       * tolScore (*double vector*)  : max flux discrepancy of the model
%                       * xopt (*cell*)               : simulated metabolite and enzyme concentrations for each experimental condition
%                       * simFluxes (*cell*)          : simulated fluxes by the model
%               * replenished particles (*int*)   : number of models generated in the current population
%
% .. Authors:
%       - Pedro Saa     2016 original code
%       - Marta Matos	2019 refactored code and added check for valid
%                       models

% 1. Load information
popIdx   = 1;
ensemble = loadEnsembleStructure(inputFile);
ensemble.eigThreshold = eigThreshold;

% This is used to always specify the same seed
%  for the random number generator when in testing in 
%  parallel mode
if nargin == 5
    ensemble.testing = testing;
else
    ensemble.testing = false;
end

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);
addKineticFxnsToPath(ensemble);


% 3. Sample fluxes and Gibbs energies
ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples);


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
tolScore = zeros(ensemble.replenishedParticles(popIdx), ensemble.numConditions);
strucIdx = zeros(ensemble.replenishedParticles(popIdx), 1);
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
            if ensemble.testing == true
                rng(1+ix);
            else
                rng(sum(clock)+ix)                                             % This is necessary to avoid generating the same results each time
            end
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
        xopt = xopt(validModelList);
        tolScore = tolScore(validModelList);
        simFluxes = simFluxes(validModelList);            
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
            tolScore(ix,:) = tolScoreTemp;
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
    clearvars -except ensemble popIdx modelID outputFile
    save(outputFile);
else
    disp('No valid models were sampled.');
end

end

