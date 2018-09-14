% Example 1. Sample a kinetic model of the MEP pathway (reference point)
%--------------------------------------------------------------------------
% Executes GRASP workflow
%
% Inputs:       (-)
%
% Outputs:      (-)
%--------------------- Pedro Saa 2017 -------------------------------------
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns');

% 1. Load information
iter     = 1;
popIdx   = 1;

ensemble = loadEnsembleStructure('input_test/MEP_example');           % Here the test case HMP pathway model is chosen

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
            rng('shuffle');
            [models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = initialSampler(ensemble);
        end
        delete(gcp);
    else
        for ix = 1:ensemble.numParticles
            [models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = initialSampler(ensemble);
        end
    end
    
    % In the ORACLE mode we are only interested in the models
else
    if ensemble.parallel
        parpool(ensemble.numCores);																								% Initiate parallel pool and run parallel foor loop
        parfor ix = 1:ensemble.numParticles
            rng('shuffle');
            models(ix) = initialSampler(ensemble);
        end
        delete(gcp);
    else
        for ix = 1:ensemble.numParticles
            models(ix) = initialSampler(ensemble);
        end
    end
end

% Append sampling results
ensemble.populations(1).strucIdx  = strucIdx;                                                                           % model structures
ensemble.populations(1).tolScore  = tolScore;                                                                           % tolerance score
ensemble.populations(1).xopt      = xopt;                                                                               % optimal value found
ensemble.populations(1).simFluxes = simFluxes;                                                                          % simulated fluxes
ensemble.populations(1).models    = models;                                                                             % model particles
clearvars -except ensemble popIdx iter

save('output_test/ensembleSMC_rejection_MEP.mat');

