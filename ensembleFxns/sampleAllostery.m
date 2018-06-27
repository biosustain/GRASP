function [ensemble, models] = sampleAllostery(ensemble, models, strucIdx)
%SAMPLEALLOSTERY Summary of this function goes here
%   Detailed explanation goes here


ensemble.reactionFluxAllosteric =  ensemble.fluxRef;

for activRxnIdx = 1:numel(ensemble.kinActRxns)  

    reactionFlux = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));
    if ensemble.allosteric{strucIdx}(ensemble.kinActRxns(activRxnIdx))

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
        ensemble.reactionFluxAllosteric{ensemble.kinActRxns(activRxnIdx)} = reactionFlux;
        
    end
    
end

end

