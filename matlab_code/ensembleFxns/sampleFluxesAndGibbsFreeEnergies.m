function ensemble = sampleFluxesAndGibbsFreeEnergies(ensemble,maxNumberOfSamples,priorType)

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
fluxPoints = generalHR(fluxAeq,fluxLB,fluxUB,fluxX0,nSamples,nSteps,nDiscard,priorType);
ensemble.fluxPoints = fluxPoints(:,end-maxNumberOfSamples+1:end);

assert(all(all(abs(ensemble.Sred * ensemble.fluxPoints) <10^-8)), "Not all sampled flux points lead to S.v = 0");

disp('Sampling Gibbs energies')

% 2) Sample thermodynamic features
RT           = 8.314*298.15/1e3;  % [kJ/mol]
[m, n]       = size(ensemble.Sthermo);
thermoLB     = [ensemble.gibbsRanges(:,1);ensemble.DGfStdRange(:,1);ensemble.lnMetRanges(:,1)];
thermoUB     = [ensemble.gibbsRanges(:,2);ensemble.DGfStdRange(:,2);ensemble.lnMetRanges(:,2)];
thermoAeq    = [eye(size(ensemble.Sthermo,2)),-ensemble.Sthermo',-RT*ensemble.Sthermo'];
thermoX0     = ensemble.initialTMFAPoint(numel(fluxX0)+1:end);
thermoPoints = generalHR(thermoAeq,thermoLB,thermoUB,thermoX0,nSamples,nSteps,nDiscard,priorType);
ensemble.gibbsEnergies = thermoPoints(1:n,end-maxNumberOfSamples+1:end);
ensemble.metConcRef = exp(thermoPoints(n+m+1:end,end-maxNumberOfSamples+1:end));

% Check that everything is consistent
Nint = null(ensemble.Sthermo,'r');

assert(all(all(abs(Nint' * ensemble.gibbsEnergies) < 1e-6)), 'The sum of dGs for reactions involved in closed loops is not zero for all.');
assert(all(all(thermoAeq * thermoPoints(:,maxNumberOfSamples+1:end) < 1e-6)), 'Not all thermodynamic dynamic points are valid.');
assert(sum(sum(sign(ensemble.fluxPoints(ensemble.idxNotExch,:)) + sign(ensemble.gibbsEnergies))) == 0, 'There seem to be dG values inconsistent with the respective flux values. Make sure that the fluxes standard deviaton is not zero in measRates.');


disp('Fluxes and Gibbs energies successfully sampled');
end
