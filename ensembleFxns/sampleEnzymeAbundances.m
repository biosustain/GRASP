function [models] = sampleEnzymeAbundances(ensemble,models,strucIdx)
%--------------------------------------------------------------------------
% Function used to sample enzyme abundances
%
%------------------------Pedro Saa 2016, Marta Matos 2018------------------
for activRxnIdx = 1:numel(ensemble.kinActRxns)
    
    if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
         strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
       
        promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
        alphaEnzymesR  = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;

        % If it's a promiscuous reaction and it's not the first one in the set of promiscuous reactions
        if size(promiscRxnsList,1) > 0 && ensemble.kinActRxns(activRxnIdx) ~= promiscRxnsList(1)                    % Get the abundances from the first reaction in the set
            randomEnzymesR = models(1).rxnParams(promiscRxnsList(1)).enzymeAbundances';
        else
            randomEnzymesR = randg(alphaEnzymesR');                                                                 % Sample from the gamma distribution with parameter alpha
            randomEnzymesR = randomEnzymesR/sum(randomEnzymesR);
        end
            models(1).rxnParams(activRxnIdx).enzymeAbundances = randomEnzymesR';
    end
end
end