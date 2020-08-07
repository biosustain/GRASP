function [ensemble, models] = sampleAllostery(ensemble,models,strucIdx)
% Function used to sample allosteric fluxes and parameters for each 
% reaction.
%
%
% USAGE:
%
%    [ensemble, models] = sampleAllostery(ensemble, models, strucIdx)
%
% INPUT:
%    ensemble (struct):       model ensemble, see buildEnsemble for fields description
%    models (struct):         model, see initialSampler for fields description
%    strucIdx (int):          number of the model structure considered
%
% OUTPUT:
%    ensemble (struct):  initialized model ensemble, see buildEnsemble for fields description
%    models (struct):    model structure with added allosteric parameters, see initialSampler for fields description
%
% .. Authors:
%       - Pedro Saa         2016 original code
%       - Marta Matos       2018 generalized it for promiscuous reactions 
%                           and random mechanisms

RT       = 8.314*298.15/1e3;                                               % gas constant times the absolute temperature (298.15 K)
DDG_max  = 0;															   % max. Gibbs free energy difference of conformation (kJ/mol)
DDG_min  = -25;															   % min. Gibbs free energy difference of conformation (kJ/mol) ~ 1 ATP hydrolysis
bind_max = 50;															   % max. fraction of Keff/Eff
bind_min = .1;															   % min. fraction of Keff/Eff

ensemble.reactionFluxAllosteric =  models(1).refFlux;

for activRxnIdx = 1:numel(ensemble.kinActRxns)
    if ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')&&...
		~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'fixedExchange')&&...
        ~strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
    
        reactionFlux =  models(1).refFlux(ensemble.kinActRxns(activRxnIdx));
        if ensemble.allosteric{strucIdx}(ensemble.kinActRxns(activRxnIdx))

            randomEnzymesR = models(1).rxnParams(activRxnIdx).enzymeAbundances';

            % I. Sample allosteric parameters
            allostericParams = mvnrnd(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).muAllostericParams,...
                ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams);
            allostericParams = exp(allostericParams)./(1 + exp(allostericParams));									    % Back-transform auxiliary parameters
            models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).allostericFactors = allostericParams;                 % allosteric factors used in the calculation of L, KposEff & KnegEff
            subunits = ensemble.subunits{strucIdx}(ensemble.kinActRxns(activRxnIdx));

            % II. Back-calculate kinetic parameters (order: L, posEff & negEff)
            DDG = DDG_min + allostericParams(1)*(DDG_max - DDG_min);
            L   = exp(-DDG/RT);
            models(1).rxnParams(activRxnIdx).L = L;                                                                     % Save allosteric constant

            % III. Compute regulatoy term
            Q = L*(randomEnzymesR(1))^subunits;
            if (numel(allostericParams)>1)																				% There are possitive and/or negative effectors
                posEffectors = ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};                       % Possitive effectors
                negEffectors = ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};                       % Negative effectors

                % IV. Back-calculate effector parameters
                Keff = bind_min*exp(allostericParams(2:end)*log(bind_max/bind_min));

                % Both effectors present
                if ~isempty(posEffectors) && ~isempty(negEffectors)
                    KposEff = Keff(1:max(size(posEffectors)));
                    KnegEff = Keff(max(size(posEffectors))+1:max(size(posEffectors))+max(size(negEffectors)));
                    models(1).rxnParams(activRxnIdx).KposEff = KposEff;							       		            % transpose params
                    models(1).rxnParams(activRxnIdx).KnegEff = KnegEff;
                    Q = Q*((1 + sum(1./KnegEff))/((1 + sum(1./KposEff))))^subunits;
                elseif ~isempty(posEffectors)
                    models(1).rxnParams(activRxnIdx).KposEff = Keff;							                        % transpose params
                    Q = Q*((1 + sum(1./Keff)))^-subunits;
                elseif ~isempty(negEffectors)
                    models(1).rxnParams(activRxnIdx).KnegEff = Keff;                        							% transpose params
                    Q = Q*((1 + sum(1./Keff)))^subunits;
                end
            end
            regContribution = 1/(1 + Q);

            % V. Calculate regulatory and catalytic contributions at the reference state
            catContribution = reactionFlux/regContribution;
            reactionFlux    = catContribution/subunits;
            ensemble.reactionFluxAllosteric(ensemble.kinActRxns(activRxnIdx)) = reactionFlux;

        end
    end
end

end
