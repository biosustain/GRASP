function [ensemble,currTol,exitFlag] = setNewTolerance(ensemble,popIdx,iter)
%--------------------------------------------------------------------------
% Compute tolerance and particles weights for the current step
%
% Inputs:       ensemble (structure) , workerIdx (double)
%
% Outputs:      initial sampled ensemble structure
%--------------------- Pedro Saa 2016 -------------------------------------
% Minimum model probability
exitFlag          = 0;   														% flag indicating status of sampler (0 = not finished, 1 = final tol reached, 2 = acceptance rate too low)
minModelStrucProb = 1e-2;                           							% Set the minimum probability for each model structure
ensemble.populations(popIdx).structureWeights = ones(1,ensemble.numStruct);

% Compute model structure weights
for ix = 1:ensemble.numStruct
	
	% Drop model structures poorly represented in the posterior
    ixStrucParticles = find(ensemble.populations(popIdx).strucIdx==ix);
	probModelIdx     = sum(ensemble.populations(popIdx).weights(ixStrucParticles));
    if (probModelIdx < minModelStrucProb)    
        ensemble.populations(popIdx).structureWeights(1,ix) = 0;
        
        % Drop particles with zero weights (based on the previous step)
        ensemble.populations(popIdx).weights(ixStrucParticles)   = [];
        ensemble.populations(popIdx).strucIdx(ixStrucParticles)  = [];
		ensemble.populations(popIdx).models(ixStrucParticles)    = [];
		ensemble.populations(popIdx).tolScore(ixStrucParticles)  = [];
		ensemble.populations(popIdx).xopt(ixStrucParticles)      = [];
		ensemble.populations(popIdx).simFluxes(ixStrucParticles) = [];
    end
end

% Re-normalize weights of the population
ensemble.populations(popIdx).weights = robustWeightNorm(ensemble.populations(popIdx).weights);

% Find the alpha-percentile of dropped particles
currTol = prctile(ensemble.populations(popIdx).tolScore,ensemble.alphaAlive);
currTol = max([currTol,ensemble.tolerance(end)]);

% Find particles above the desired tolScore
droppedParticles = find(ensemble.populations(popIdx).tolScore>currTol);

% Remove particles based on the tolerance
ensemble.populations(popIdx).weights(droppedParticles)   = [];
ensemble.populations(popIdx).strucIdx(droppedParticles)  = [];
ensemble.populations(popIdx).models(droppedParticles)    = [];
ensemble.populations(popIdx).tolScore(droppedParticles)  = [];
ensemble.populations(popIdx).xopt(droppedParticles)      = [];
ensemble.populations(popIdx).simFluxes(droppedParticles) = [];

% Re-normalize weights of the population
ensemble.populations(popIdx).weights = robustWeightNorm(ensemble.populations(popIdx).weights);

% 1. Extract valid particles from the previous population
for idxStruc = 1:ensemble.numStruct       
    
    % Find particles belonging to this model structure
    ixStrucParticles = find(ensemble.populations(popIdx).strucIdx==idxStruc);
    
	% Check that the model structure is valid
    if isempty(ixStrucParticles)
		ensemble.populations(popIdx).structureWeights(1,idxStruc) = 0;
		continue; 
	end 
	
	% Drop model structures poorly represented in the posterior
	probModelIdx = sum(ensemble.populations(popIdx).weights(ixStrucParticles));
    if (probModelIdx < minModelStrucProb)    
        ensemble.populations(popIdx).structureWeights(1,ix) = 0;
        
        % Drop particles with zero weights (based on the previous step)
        ensemble.populations(popIdx).weights(ixStrucParticles)   = [];
        ensemble.populations(popIdx).strucIdx(ixStrucParticles)  = [];
		ensemble.populations(popIdx).models(ixStrucParticles)    = [];
		ensemble.populations(popIdx).tolScore(ixStrucParticles)  = [];
		ensemble.populations(popIdx).xopt(ixStrucParticles)      = [];
		ensemble.populations(popIdx).simFluxes(ixStrucParticles) = [];
    end
end

% Sort tolerance accordingly
if (currTol~=ensemble.tolerance(end))
    ensemble.tolerance = [ensemble.tolerance(1:end-1),currTol,ensemble.tolerance(end)];
else
	exitFlag = 1;											
end

% Normalize particle and model weights
ensemble.populations(popIdx).weights = robustWeightNorm(ensemble.populations(popIdx).weights);											% Normalize particle weights
ensemble.ESS(iter-1)                 = 1/sum(ensemble.populations(popIdx).weights.^2);													% Compute effective sample size
ensemble.aliveParticles(iter-1)      = sum(ensemble.populations(popIdx).weights~=0);													% Determine number of alive particles from the previous population

% Determine adaptive weights based on robust model selection algorithm
if ensemble.robustModSel
	for ix = 1:ensemble.numStruct		
		ensemble.populations(popIdx).structureWeights(ix) = sum(ensemble.populations(popIdx).weights(ensemble.populations(popIdx).strucIdx==ix));
	end
end	

% Display parameters for the next population
disp(['Effective tolerance: ',num2str(currTol),'. Model weights for next step: ',mat2str(ensemble.populations(popIdx).structureWeights)]);
disp(['Number of alive particles for resampling: ',num2str(ensemble.aliveParticles(iter-1)),'.']);
disp(['Effective Sample Size for this iteration: ',num2str(ensemble.ESS(iter-1)),'.']);