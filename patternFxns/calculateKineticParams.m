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
elemFluxVector = Nelem*branchFactor / max(Nelem*branchFactor);

%elemFluxVector = Nelem*branchFactor/sum(Nelem*branchFactor);
%disp('Nelem');
%disp(Nelem);
%disp('branchFactor');
%disp(branchFactor);
%disp('Nelem*branchFactor');
%disp(Nelem*branchFactor);
%disp('max(Nelem*branchFactor)');
%disp(max(Nelem*branchFactor));
%disp('elemFluxVector');
%disp(elemFluxVector);


% For promiscuous reactions, set inhib entries to nan
if size(Nelem,2) > 1 && sum(sum(Nelem)) <= size(Nelem,1)

    inhibEntries = find(all(revCalIrrev==1, 2));
    if size(inhibEntries) > 0
        nTracks = size(revCalIrrev, 2);
        
        for j = 1:nTracks;
            for i = 1:size(inhibEntries)
                inhibEntry = inhibEntries(i);
                if revCalIrrev(inhibEntry-1, j) ~= 1 && revCalIrrev(inhibEntry+1, j) ~=1
                    revCalIrrev(inhibEntry, j) = nan;
                end
            end
        end
    end
         
end

% If the proposed branch vector is OK continue
revCal     = [(revCalIrrev.*elemFluxVector)';(revTemp.*revCalIrrev.*elemFluxVector)'];
elemenFlux = revCal(:);

%disp('elemenFlux');
%disp(elemenFlux);

%disp('modifierElemFlux');
%disp(modifierElemFlux);



% If the pattern contains modifiers (e.g. inhibitors/activators) compute
% elementary fluxes based on the elem flux fraction
if any(isnan(elemenFlux))
    elemenFlux(isnan(elemenFlux)) = modifierElemFlux;
end


%disp('reactionFlux');
%disp(reactionFlux);

%disp('elemenFlux');
%disp(elemenFlux);

%disp('enzymeVect');
%disp(enzymeVect);

% 3. Output the kinetic parameters
K = reactionFlux*elemenFlux.*(enzymeVect(:).^(-1));
%disp('K');
%disp(K);
%disp('----');