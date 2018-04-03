function ensemble = collectParallelEnsembles(ensembles,popIdx)
%--------------------------------------------------------------------------
% Collect parallel ensembles
%
% Inputs:       ensembles (cell array of structures)
%
% Outputs:      ensemble structure
%--------------------- Pedro Saa 2016 -------------------------------------
ensemble       = ensembles{1};                                             % initial assignation is trivial
currTime       = ensemble.time(popIdx);
currAcceptance = ensemble.acceptanceRate(popIdx);

% Loop through the remaining structures
if numel(ensembles)>1
    for ix = 2:numel(ensembles)
        
        % Assign temporal structure
        ensembleTemp = ensembles{ix};
        
        % Add new particles to the existing structure
        ensemble.populations(popIdx).strucIdx  = [ensemble.populations(popIdx).strucIdx;ensembleTemp.populations(popIdx).strucIdx];          % model structures
        ensemble.populations(popIdx).tolScore  = [ensemble.populations(popIdx).tolScore;ensembleTemp.populations(popIdx).tolScore];          % tolerance score        
        ensemble.populations(popIdx).xopt      = [ensemble.populations(popIdx).xopt;ensembleTemp.populations(popIdx).xopt];                  % optimal value found
        ensemble.populations(popIdx).simFluxes = [ensemble.populations(popIdx).simFluxes;ensembleTemp.populations(popIdx).simFluxes];        % simulated fluxes
        ensemble.populations(popIdx).models    = [ensemble.populations(popIdx).models,ensembleTemp.populations(popIdx).models];              % model particles
        ensemble.populations(popIdx).weights   = [ensemble.populations(popIdx).weights;ensembleTemp.populations(popIdx).weights];            % importance weights (these are all zero)
        
        % Assign new values for time and acceptance
        currTime       = [currTime,ensembleTemp.time(end)];
        currAcceptance = [currAcceptance,ensembleTemp.acceptanceRate(end)];        
    end
end

% Extract models from previous population (from the second population onwards)
if (popIdx>1)
    activeParticles                        = find(ensemble.populations(popIdx-1).weights~=0);
    disp(['Number of active particles in the prev population ',num2str(numel(activeParticles))])
    ensemble.populations(popIdx).strucIdx  = [ensemble.populations(popIdx).strucIdx;ensemble.populations(popIdx-1).strucIdx(activeParticles)];        % model structures
    ensemble.populations(popIdx).tolScore  = [ensemble.populations(popIdx).tolScore;ensemble.populations(popIdx-1).tolScore(activeParticles)];        % tolerance score
    ensemble.populations(popIdx).xopt      = [ensemble.populations(popIdx).xopt;ensemble.populations(popIdx-1).xopt(activeParticles)];                % optimal value found
    ensemble.populations(popIdx).simFluxes = [ensemble.populations(popIdx).simFluxes;ensemble.populations(popIdx-1).simFluxes(activeParticles)];      % simulated fluxes
    ensemble.populations(popIdx).models    = [ensemble.populations(popIdx).models,ensemble.populations(popIdx-1).models(activeParticles)];            % model particles
    ensemble.populations(popIdx).weights   = [ensemble.populations(popIdx).weights;ensemble.populations(popIdx-1).weights(activeParticles)];          % importance weights (these are non-zero)
end

% Re-compute statistics report
ensemble.time(popIdx)           = mean(currTime);
ensemble.acceptanceRate(popIdx) = mean(currAcceptance);
ensemble.numParticles           = numel(ensemble.populations(popIdx).weights);                                                                           % Update number of particles