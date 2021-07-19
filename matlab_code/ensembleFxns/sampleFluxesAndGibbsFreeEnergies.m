function ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples)

% Define sampling parameters and sample
nDiscard = 1e5;
nSamples = nDiscard+maxNumberOfSamples;
nSteps   = 1e2;

disp('Sampling fluxes')
% 1) Sample fluxes indenpendently as the problem is uncoupled
fluxLB     = ensemble.fluxRanges(:,1);
fluxUB     = ensemble.fluxRanges(:,2);
fluxX0     = ensemble.initialTMFAPoint(1:size(fluxLB,1));
fluxAeq    = ensemble.Sflux;
fluxPrior  = ensemble.fluxPrior;

if ~ensemble.parallel
    fluxPoints = generalHR(fluxAeq,fluxLB,fluxUB,fluxX0,nSamples,nSteps,nDiscard,fluxPrior);
else
    fluxPoints{ensemble.numCores} = [];
    nSamplesPerCore = round(nSamples/ensemble.numCores);
    parpool(ensemble.numCores);                 % Run one Markov chain per core
    parfor ix = 1:ensemble.numCores
        if ensemble.testing == true
            rng(1+ix)                   % The results need to be reproducible when testing
        else
            rng(sum(clock)+ix)     % This is necessary to avoid generating the same results by the workers 
        end
        fluxPoints{ix} = generalHR(fluxAeq,fluxLB,fluxUB,fluxX0,nSamplesPerCore,nSteps,nDiscard,fluxPrior);
    end
    delete(gcp('nocreate'));
    fluxPoints = cell2mat(fluxPoints);
end

ensemble.fluxPoints = fluxPoints(:,end-maxNumberOfSamples+1:end);

assert(all(all(abs(ensemble.Sred * ensemble.fluxPoints) <10^-8)), "Not all sampled flux points lead to S.v = 0");


disp('Sampling Gibbs energies')

% 2) Sample thermodynamic features
RT           = 8.314*298.15/1e3;  % [kJ/mol]
[m, n]       = size(ensemble.Sthermo);
thermoLB     = [ensemble.gibbsRanges(ensemble.idxNotExch,1);ensemble.DGrStdRange(ensemble.idxNotExch,1);ensemble.lnMetRanges(:,1)];
thermoUB     = [ensemble.gibbsRanges(ensemble.idxNotExch,2);ensemble.DGrStdRange(ensemble.idxNotExch,2);ensemble.lnMetRanges(:,2)];
thermoAeq    = [eye(size(ensemble.Sthermo,2)),-eye(size(ensemble.Sthermo,2)),-RT*ensemble.Sthermo'];
thermoX0     = ensemble.initialTMFAPoint(numel(fluxX0)+1:end);
thermoPrior  = ensemble.thermoPrior;

if ~ensemble.parallel
    thermoPoints = generalHR(thermoAeq,thermoLB,thermoUB,thermoX0,nSamples,nSteps,nDiscard,thermoPrior);
else
    thermoPoints{ensemble.numCores} = [];
    parpool(ensemble.numCores);                 % Run one Markov chain per core
    parfor ix = 1:ensemble.numCores
        if ensemble.testing == true
            rng(1+ix)                   % The results need to be reproducible when testing
        else
            rng(sum(clock)+ix)     % This is necessary to avoid generating the same results by the workers 
        end                      % This is necessary to avoid generating the same results by the workers 
        thermoPoints{ix} = generalHR(thermoAeq,thermoLB,thermoUB,thermoX0,nSamplesPerCore,nSteps,nDiscard,fluxPrior);
    end
    delete(gcp('nocreate'));
    thermoPoints = cell2mat(thermoPoints);
end
ensemble.gibbsEnergies = thermoPoints(1:n,end-maxNumberOfSamples+1:end);
ensemble.metConcRef = exp(thermoPoints(n+m+1:end,end-maxNumberOfSamples+1:end));


% Check that everything is consistent
Nint = null(ensemble.Sthermo,'r');

assert(all(all(abs(Nint' * ensemble.gibbsEnergies) < 1e-4)), 'The sum of dGs for reactions involved in closed loops is not zero for all.');
assert(all(all(thermoAeq * thermoPoints(:,maxNumberOfSamples+1:end) < 1e-4)), 'Not all thermodynamic dynamic points are valid.');
assert(sum(sum(sign(ensemble.fluxPoints(ensemble.idxNotExch,:)) + sign(ensemble.gibbsEnergies))) == 0, 'There seem to be dG values inconsistent with the respective flux values. Make sure that the fluxes standard deviaton is not zero in measRates.');


disp('Fluxes and Gibbs energies successfully sampled');
end

