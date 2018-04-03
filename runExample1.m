% Example 1. Sample a kinetic model of the MEP pathway (reference point)
%--------------------------------------------------------------------------
% Executes GRASP workflow
%
% Inputs:       (-)
%
% Outputs:      (-)
%--------------------- Pedro Saa 2017 -------------------------------------
profile on
rng('default');																											% for reproducibility
delete(gcp('nocreate'));       				            																% check first that no other process is running
addpath('./patternFxns','./ensembleFxns','./reactions_MEP');

% 1. Load information
iter     = 1;
popIdx   = 1;
%ensemble = loadEnsembleStructure('models/input/melatonin_pathway_no_transporters_HMP1489');		% This line loads the model in the Excel file
ensemble = loadEnsembleStructure('models/input/MEP_example');		% This line loads the model in the Excel file
%load('ensembleStarter.mat','ensemble');

% 2. Initialize and perform rejection sampling
ensemble = initializeEnsemble(ensemble,popIdx,1);

% Check whether the job is ran in parallel
disp('Running rejection sampler. Population 1.');

% Setup folder with temp files
try
	rmdir('tempRejection','s');
	mkdir('tempRejection');
catch
	mkdir('tempRejection');
end
	
% Preallocate memory for the remaing fields in the ensemble structure
tolScore = zeros(ensemble.replenishedParticles(popIdx),1);
strucIdx = zeros(ensemble.replenishedParticles(popIdx),1);    
xopt{ensemble.replenishedParticles(popIdx),1}      = [];
simFluxes{ensemble.replenishedParticles(popIdx),1} = [];
    
% Initiate progress report	
progress = zeros(5,1);	
save progress.txt -ascii progress;    
if ensemble.parallel
		parpool(ensemble.numCores);																								% Initiate parallel pool and run parallel foor loop    
	parfor ix = 1:ensemble.numParticles
		rng('shuffle');
		[models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = rejectionSampler(ensemble);
	end
	delete(gcp);
else
	for ix = 1:ensemble.numParticles
		[models(ix),strucIdx(ix),xopt{ix},tolScore(ix),simFluxes{ix}] = rejectionSampler(ensemble);
	end
end

% Append active particles from the previous population to the new ones
ensemble.populations(1).strucIdx  = strucIdx;                                                                           % model structures
ensemble.populations(1).tolScore  = tolScore;                                                                           % tolerance score
ensemble.populations(1).xopt      = xopt;                                                                               % optimal value found
ensemble.populations(1).simFluxes = simFluxes;                                                                          % simulated fluxes
ensemble.populations(1).models    = models;                                                                             % model particles    
clearvars -except ensemble popIdx iter
save('models/output/ensembleSMC_rejection_MEP.mat');
profsave