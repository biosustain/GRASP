function K = calculateKineticParams(randomRev,forwardFlux,reactionFlux,randomEnzymes,Nelem,branchFactor,modifierElemFlux)
%--------------------------------------------------------------------------
% Function used to calculate the kinetic parameters of the reaction
%
% Inputs:    (randomRev)  feasible random reversibilities in the pattern
%          (forwardFlux)  list with the edges of the forward reactions
%         (reactionFlux)  overall reaction flux
%        (randomEnzymes)  random enzyme complexes abundances
%                (Nelem)  null basis Selem
%         (branchFactor)  branching factor ~distributes Dir(1)
%   (modifierFluxFactor)  elementary flux entries for mechanisms with modifiers
%
% Outputs:           (K)  kinetic parameters for the reaction pattern
%-----------------------Pedro Saa 2016-------------------------------------
% 1. Sort enzyme abundances
enzymeVect = [randomEnzymes(forwardFlux(:,1))';randomEnzymes(forwardFlux(:,2))'];

% 2. Calculate elementary fluxes
revTemp      = randomRev.^(sign(reactionFlux));
revCalIrrev  = (1-revTemp).^(-1);

% Compute branching flux structure
elemFluxVector = Nelem*branchFactor/max(Nelem*branchFactor);

% If the proposed branch vector is OK continue
revCal     = [(revCalIrrev.*elemFluxVector)';(revTemp.*revCalIrrev.*elemFluxVector)'];
elemenFlux = revCal(:);

% If the pattern contains modifiers (e.g. inhibitors/activators) compute
% elementary fluxes based on the elem flux fraction
if any(isnan(elemenFlux))
    elemenFlux(isnan(elemenFlux)) = modifierElemFlux;
end

% 3. Output the kinetic parameters
K = reactionFlux*elemenFlux.*(enzymeVect(:).^(-1));