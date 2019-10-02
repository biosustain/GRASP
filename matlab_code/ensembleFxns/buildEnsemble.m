function ensemble = buildEnsemble(inputFile,outputFile,maxNumberOfSamples,eigThreshold)
% Samples a kinetic model ensemble which contains only valid models.
%
% Valid models are models where:
%  - for all reactions the fluxes and respective Gibbs energies are 
%    compatible;
%  - the real part of the jacobian eigenvalues is lower than the defined 
%    threshold;
%  - the difference between the predicted flux and the reference flux
%    is negligible.
%
% It will sample models until there are n valid models, where n is
% specified in the input excel file as the number of particles.
% To avoid sampling forever, maxNumberOfSamples is defined and no more
% than maxNumberOfSamples models will be sampled, regardless of the number
% of valid models sampled.
% 
% Models with a real part of the jacobian eigenvalues greather than 
% eigThreshold are discarded.
%
% If the output directory doesn't exist it is created.
%
% USAGE:
%
%    ensemble = buildEnsemble(inputFile, outputFile, maxNumberOfSamples,
%                             eigThreshold)
%
% INPUTS:
%    inputFile (`char`):            path to input file
%    outputFile (`char`):           path to output file
%    maxNumberOfSamples (`int`):	maximum number of models to be sampled
%    eigThreshold (`double`):       threshold for positive eigenvalues' 
%                                   real part
%
% OUTPUT:
%    ensemble (`struct`):    model ensemble
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


% Check whether the job is ran in parallel
disp('Running initial sampler.');

% Setup folder with temp files
try
    rmdir('temp','s');
    mkdir('temp');
catch
    mkdir('temp');
end

% Preallocate memory for the remaing fields in the ensemble structure
tolScore = zeros(ensemble.replenishedParticles(popIdx),1);
strucIdx = zeros(ensemble.replenishedParticles(popIdx),1);
xopt{ensemble.replenishedParticles(popIdx),1}      = [];
simFluxes{ensemble.replenishedParticles(popIdx),1} = [];

% Figure out the sampling mode
if ~strcmpi(ensemble.sampler,'ORACLE')
    
    % Initiate progress report for rejection sampling
    progress = zeros(5,1);
    save progress.txt -ascii progress;
    if ensemble.parallel
        parpool(ensemble.numCores);																								% Initiate parallel pool and run parallel foor loop
        parfor ix = 1:ensemble.numParticles
            [isModelValid,models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = initialSampler(ensemble);
        end
        delete(gcp);
    else
        for ix = 1:ensemble.numParticles
            [isModelValid,models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = initialSampler(ensemble);
        end
    end
    
% In the ORACLE mode we are only interested in the models
else

    if ensemble.parallel
        sampleCount = 0;
        nValidModels = 0;
        
        while nValidModels < ensemble.numParticles
            
            if sampleCount == 0                                             % if no models have been sampled yet, sample ensemble.numParticles models
                nSamples = ensemble.numParticles;
            else                                                            % else check how many more should be sampled based on the percentage of valid models
                nSamples = 1 / (nValidModels / sampleCount)*(sampleCount - nValidModels);
                nSamples = round(nSamples);
                
                if nSamples > maxNumberOfSamples - sampleCount
                    nSamples = maxNumberOfSamples - sampleCount;
                end
            end
            
            parpool(ensemble.numCores);
            parfor ix = (sampleCount+1):(sampleCount+nSamples)
                [validModelList(ix),models(ix)] =  initialSampler(ensemble);
            end
            delete(gcp);
            
            sampleCount = sampleCount + nSamples;           
            nValidModels = size(models(validModelList), 2);
        
        end
        models = models(validModelList);
        
        
    else
        sampleCount = 1;
        ix = 1;
        
        while ix <= ensemble.numParticles && sampleCount <= maxNumberOfSamples
            
            [isModelValid,model] =  initialSampler(ensemble);
            if isModelValid
                models(ix) = model;
                ix = ix + 1;
            end
            
            sampleCount = sampleCount + 1;
        end
        
        nValidModels = ix - 1;
        sampleCount = sampleCount -1;
    end
    
    disp([newline, newline, '*** Out of ',num2str(sampleCount), ' models, ', num2str(nValidModels), ' were valid ***', newline, newline]);
end


if nValidModels > 1
    % Append sampling results
    ensemble.populations(1).strucIdx  = strucIdx;                                                                           % model structures
    ensemble.populations(1).tolScore  = tolScore;                                                                           % tolerance score
    ensemble.populations(1).xopt      = xopt;                                                                               % optimal value found
    ensemble.populations(1).simFluxes = simFluxes;                                                                          % simulated fluxes
    ensemble.populations(1).models    = models;                                                                             % model particles
    clearvars -except ensemble popIdx modelID outputFile
    save(outputFile);
end

end

