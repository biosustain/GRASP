function [ensemble] = sampleGibbsReactionEnergies(ensemble, gibbsEnergy, strucIdx)
%--------------------------------------------------------------------------
% Function used to sample Gibbs energies for each reaction
%
%------------------------Pedro Saa 2016, Marta Matos 2018------------------


ensemble.gibbsTemp = cell(size(ensemble.thermoActive));
thermoCounter = 1;
for activRxnIdx = 1:numel(ensemble.kinActRxns)        

    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))

        % Determine gibbs free energy of reaction
        if ismember(ensemble.kinActRxns(activRxnIdx),ensemble.thermoActive)
            ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)}     = gibbsEnergy(thermoCounter);
            thermoCounter = thermoCounter+1;
        elseif (ensemble.gibbsRanges(activRxnIdx,1)==ensemble.gibbsRanges(activRxnIdx,2))
            ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)}     = ensemble.gibbsRanges(activRxnIdx,1);
        end

    end
end

end

