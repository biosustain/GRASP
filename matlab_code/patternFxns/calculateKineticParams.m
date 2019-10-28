function K = calculateKineticParams(randomRev,forwardFlux,reactionFlux,randomEnzymes,Nelem,branchFactor,modifierElemFlux,rxnIsPromiscuous)
% Function used to calculate the kinetic parameters of the reaction.
%
%
% USAGE:
%
%    K = calculateKineticParams(randomRev, forwardFlux, reactionFlux, randomEnzymes, Nelem, branchFactor, modifierElemFlux, rxnIsPromiscuous)
%
% INPUT:
%    randomRev (double vector):           feasible random reversibilities in  the pattern 
%    forwardFlux (int matrix):            list with the edges of the forward reactions
%    reactionFlux (double vector):        overall reaction flux
%    randomEnzymes (double vector):       random enzyme complexes abundances
%    Nelem (int matrix):                  null basis Selem
%    branchFactor (int vector):           branching factor ~distributes  Dir(1) 
%    modifierElemFlux (double vector):	  elementary flux entries for mechanisms with modifiers
%    rxnIsPromiscuous (logical):          boolean specifying whether or not reaction is promiscuous
%
% OUTPUT:
%    K  (double vector):	kinetic parameters for the reaction pattern
%
% .. Authors:
%       - Pedro Saa    2016 original code
%       - Marta Matos  2019 extended for promiscuous enzymes and for
%                      reverse reactions



% 1. Sort enzyme abundances
enzymeVect = [randomEnzymes(forwardFlux(:,1))';randomEnzymes(forwardFlux(:,2))'];

% 2. Calculate elementary fluxes
if rxnIsPromiscuous
    revTemp      = sum(Nelem.*randomRev.^(sign(reactionFlux*branchFactor)), 2);
    revTemp(find(revTemp == 0)) = 1;
else
    revTemp      = randomRev.^(sign(reactionFlux));
end

revCalIrrev  = (1-revTemp).^(-1);


% Compute branching flux structure
elemFluxVector = Nelem*branchFactor';


% For promiscuous reactions, set inhib entries to nan
if size(Nelem,2) > 1 && sum(sum(Nelem)) <= size(Nelem,1)

    inhibEntries = find(all(revCalIrrev==1, 2));
    if size(inhibEntries) > 0
        nTracks = size(revCalIrrev, 2);
        
        for j = 1:nTracks
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


% If the pattern contains modifiers (e.g. inhibitors/activators) compute
% elementary fluxes based on the elem flux fraction
if any(isnan(elemenFlux))
    elemenFlux(isnan(elemenFlux)) = sign(reactionFlux)*modifierElemFlux;
end

  
% 3. Output the kinetic parameters
K = reactionFlux*elemenFlux.*(enzymeVect(:).^(-1));

assert(all(K >= 0), 'There are negative kinetic parameters, good luck! :)');
