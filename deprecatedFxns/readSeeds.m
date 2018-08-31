function readSeeds(folderName,fileName,flag)
rng('default');                                         					% for reproducibility
addpath('ensembleFxns')
load(fileName);
method = 1;

% From rejection samples
if (flag==0)	
	iter   = 1;
	popIdx = 1;
	files  = dir(folderName);
	ensemble.numParticles = numel(files)-2;
	ensemble = initializeEnsemble(ensemble,popIdx,1);
	ensemble.alphaAlive = 30;

	% Loop for reading files
	counter = 3;
	filesLoaded = 0;
	while (counter<=numel(files))
		try
			load([folderName,'/',files(counter).name]);			
			strucIdxPop(counter-2,1)  = strucIdx;
			tolScorePop(counter-2,1)  = tolScore;		
			xoptPop{counter-2,1}      = xopt;
			simFluxesPop{counter-2,1} = simulatedFlux;
			if (counter==3)
				modelsPop = models;
			else
				modelsPop = [modelsPop,models];
			end		
			counter = counter+1;
			filesLoaded = filesLoaded + 1;
		catch			
			disp(['Corrupted particle id: ',files(counter).name]);
			counter = counter+1;
		end
	end
	disp(['Number of particles loaded: ',num2str(filesLoaded)]);		
	if (filesLoaded~=ensemble.numParticles)
		disp('Remove corrupted particles and re-run this script.');
		return;
	end
	
	% Assign data
	ensemble.populations(1).strucIdx  = strucIdxPop;        					% model structures
	ensemble.populations(1).tolScore  = tolScorePop;        					% tolerance score
	ensemble.populations(1).xopt      = xoptPop;                				% optimal value found
	ensemble.populations(1).simFluxes = simFluxesPop;      						% simulated fluxes
	ensemble.populations(1).models    = modelsPop;   		 					% model particles
	ensemble.numParticles             = numel(ensemble.populations(1).strucIdx);
	ensemble.populations(1).weights   = ones(ensemble.numParticles,1)/ensemble.numParticles;
	clearvars -except ensemble popIdx method iter
	save('ensembleSMC_rejection.mat');

	% Sequential sampler
else	
	
	% General iteration counter
	iter      = iter+1;    
	files     = dir(folderName);
	threshold = 99;
	ensemble.alphaAlive = 75;
	ensemble.numParticles = numel(files)-2;	
		
	% Prepare population for next step		
    [ensemble,currTol,exitFlag] = setNewTolerance(ensemble,popIdx,iter);                                                     % compute weights and tolerance

	% Update population and Kernel parameters
	popIdx   = popIdx+1;
	ensemble = computeKernelParams(ensemble,popIdx,method);                                                                  % compute parameters of the kernel

	% Compute weights for the active particles	
	prevParams = getPopulationParams(ensemble,popIdx);
    log_particleWeights = zeros(ensemble.aliveParticles(iter-1),1);
	parpool(ensemble.numCores);
	parfor ix = 1:ensemble.aliveParticles(iter-1)
		model    = ensemble.populations(popIdx-1).models(ix);
		strucIdx = ensemble.populations(popIdx-1).strucIdx(ix);			
		log_particleWeights(ix) = computeParticleWeightsVectorized(model,strucIdx,ensemble,popIdx,prevParams);
	end		
	delete(gcp);
	
	% Determine statistics for the log_weights
	log_weightTol = prctile(log_particleWeights,threshold);
	disp('Un-normalized weights updated.')
	
	% Accept/reject particles below the PRC threshold
	accept_prob = exp(min([zeros(size(log_particleWeights)),log_particleWeights - log_weightTol],1));
	accept_prob = accept_prob(:,2);
	activeParticles = (rand(size(log_particleWeights)) < accept_prob);
	log_particleWeights(activeParticles)  = log_particleWeights(activeParticles) - log(accept_prob(activeParticles));
	log_particleWeights(~activeParticles) = [];
	
	% Determine particles with high weights
	ensemble.populations(popIdx).strucIdx  = ensemble.populations(popIdx-1).strucIdx(activeParticles);       % model structures
	ensemble.populations(popIdx).tolScore  = ensemble.populations(popIdx-1).tolScore(activeParticles);       % tolerance score
	ensemble.populations(popIdx).xopt      = ensemble.populations(popIdx-1).xopt(activeParticles);           % optimal value found
	ensemble.populations(popIdx).simFluxes = ensemble.populations(popIdx-1).simFluxes(activeParticles);      % simulated fluxes
	ensemble.populations(popIdx).models    = ensemble.populations(popIdx-1).models(activeParticles);   		 % model particles
	ensemble.populations(popIdx).weights   = log_particleWeights;          				 					 % set importance weights of the alive particles
	% ensemble.replenishedParticles(popIdx)  = ensemble.numParticles-numel(activeParticles);  			     % Initialize number of particles, i.e., particles to be replenished = (#total-#NoWeights)/#cores

	% Load previous particles
	filesLoaded   = 0;
	counterFailed = 0;
	for counter = 3:numel(files)
		try
			load([folderName,'/',files(counter).name]);			
			filesLoaded = filesLoaded + 1;
			log_particleWeightPop(filesLoaded,1) = log_particleWeight;
			strucIdxPop(filesLoaded,1)           = strucIdx;
			tolScorePop(filesLoaded,1)           = tolScore;
			xoptPop{filesLoaded,1}               = xopt;
			simFluxesPop{filesLoaded,1}          = simulatedFlux;						
			if (filesLoaded==1)
				modelsPop = models;
			else
				modelsPop = [modelsPop,models];
			end			
		catch
			disp(['Corrupted particle id: ',files(counter).name]);
			counterFailed = counterFailed + 1;
		end
	end

	disp(['Number of particles loaded: ',num2str(filesLoaded)]);		
	if (counterFailed>0)
		disp('Remove corrupted particles and re-run this script.');
		return;
	end

	% Append new particles				
	ensemble.populations(popIdx).strucIdx  = [ensemble.populations(popIdx).strucIdx(:);strucIdxPop];   				   % model structures
	ensemble.populations(popIdx).tolScore  = [ensemble.populations(popIdx).tolScore(:);tolScorePop];                   % tolerance score
	ensemble.populations(popIdx).xopt      = [ensemble.populations(popIdx).xopt(:);xoptPop];                           % optimal value found
	ensemble.populations(popIdx).simFluxes = [ensemble.populations(popIdx).simFluxes(:);simFluxesPop];                 % simulated fluxes
	ensemble.populations(popIdx).models    = [ensemble.populations(popIdx).models,modelsPop]; 		                   % model particles
	ensemble.populations(popIdx).weights   = [ensemble.populations(popIdx).weights(:);log_particleWeightPop];          % log importance weights
	ensemble.replenishedParticles(popIdx)  = filesLoaded;
	
	% Re-normalize the un-normalized logarithmic weights
	ensemble.populations(popIdx).weights = ensemble.populations(popIdx).weights - max(ensemble.populations(popIdx).weights);
	ensemble.populations(popIdx).weights = exp(ensemble.populations(popIdx).weights)/sum(exp(ensemble.populations(popIdx).weights));
	ensemble.numParticles	             = numel(ensemble.populations(popIdx).weights);
	
	% Display ESS after loading the particles
	disp(['ESS after this iteration: ',num2str(1/sum(ensemble.populations(popIdx).weights.^2))]);
	
	% Perform population swap
    if (iter>2)
        ensemble.populations(2) = [];                
        popIdx = 2;
    end            
    clearvars -except ensemble popIdx exitFlag method iter currTol 														% cleanup workspace
	save(['ensembleSMC_',num2str(iter)]);                                                                 				% save results
end