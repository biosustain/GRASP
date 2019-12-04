function calculateMMCurves(outputFolder, ensemble, numModels, structIdx, saturatingConc, substrateRange, rxnList)
% Calculates Michaelis-Menten curves for enzymatic reactions in the model
% ensemble (i.e., reactions that are not massAction, diffusion,
% freeExchange, and fixedExchange).
%
% This basically implies setting product concentrations to 0, cosubstrate
% concentrations to a saturating value (so that they are not limiting), and
% varying the concentration of the substrate of interest from a very low 
% concentration value to a saturating one. 
%
% The results are written in files that go into outputFolder.
%
% The user can also specify the saturating concentration to be considered
% as well as the range of concentrations for which the flux 
%
% The user can further specify, for which reactions will the MM curves be
% calculated.
%
% USAGE:
%
%    mcaResults = controlAndResponseAnalysis(ensemble, saveResMatrices, strucIdx)
%
% INPUT:
%    outputFolder (char):         path to folder where output files will be written
%    ensemble (struct):           model ensemble, see buildEnsemble for fields description
%    numModels (int):             number of models in the ensemble to consider
%
% OPTIONAL INPUT:
%    structIdx (int):             structure ID, default: 1
%    saturatingConc (int):        saturation concentration (that would lead to Vmax), default: 10 mol/L
%    substrateRange (vector):     vector of concentrations for which the flux will be calculated, default: [10^-9, 10] mol/L
%    rxnList (vector):            reactions for which Vmax will be calculated
%
% OUTPUT:
%    None
%
% .. Authors:
%       - Marta Matos   2019 original code

if (nargin < 3)
    error('At least 3 arguments need to be specified: outputFolder, ensemble, and numModels');
end

if (nargin < 4)
	structIdx = 1;
end

if (nargin < 5) || isempty(saturatingConc)
    saturatingConc = 10^4; % in mmol/L
end

if (nargin < 6) || isempty(substrateRange)
    substrateRange = logspace(-9, 4);  % in mmol/L
end

if (nargin < 7)
	rxnList = 1:numel(ensemble.rxns);
end

if numModels > numel(ensemble.populations.models)
    numModels = numel(ensemble.populations.models);
end


for rxnI = rxnList
    disp(['Current reaction: ', num2str(rxnI)]);
    
    if sum(~ismember({'massAction', 'fixedExchange', 'freeExchange', 'diffusion'}, ensemble.rxnMechanisms{structIdx}{rxnI})) == 4 
        
        % if the reaction is allosteric we are only interested on the
        % catalytic part.
        if exist([ensemble.rxns{rxnI}, num2str(structIdx), 'Catalytic']) == 2
            rateLawFxn = str2func([ensemble.rxns{rxnI}, num2str(structIdx), 'Catalytic']);
        else
            rateLawFxn = str2func([ensemble.rxns{rxnI}, num2str(structIdx)]);
        end

        stoicSubsInd = find(ensemble.S(:, rxnI) < 0);        
        nSubs = numel(ensemble.subOrder{structIdx}{rxnI});
        nProds = numel(ensemble.prodOrder{structIdx}{rxnI});
        nInhib = numel(ensemble.inhibitors{structIdx}{rxnI});
        nActiv = numel(ensemble.activators{structIdx}{rxnI});

        for j=1:numel(stoicSubsInd)
            subI = stoicSubsInd(j);
            
            subList = [];
            vList = [];
            modelList = [];
            
            % catch: concentrations of substrates for promiscuous
            % reactions must be zero instead of saturating.
            subsConc = zeros(nSubs, 1);
            coSubsInd = find(ismember(ensemble.subOrder{structIdx}{rxnI}, ensemble.mets(stoicSubsInd)));
            
            % met position in ensemble.subOrder
            subOrderPos = find(ismember(ensemble.subOrder{structIdx}{rxnI}, ensemble.mets(subI)));

            % in case a metabolite is part of the reaction but not included
            % in the mechanism
            if isempty(subOrderPos)
                continue
            end

            % ensemble.subOrder met position in ensemble.mets
            subOrderInd = [];
            for entry=1:numel(ensemble.subOrder{structIdx}{rxnI})
                ind = find(ismember(ensemble.mets, ensemble.subOrder{structIdx}{rxnI}{entry}));
                if ismember(ind, stoicSubsInd)
                    subOrderInd = [subOrderInd, ind];
                end
            end
            
            for modelI=1:numModels    

                K = ensemble.populations.models(modelI).rxnParams(rxnI).kineticParams;
                allRefConcs = ensemble.populations.models(modelI).metConcRef(subOrderInd) * 10^3;
                                
                subsConc(coSubsInd) = saturatingConc ./ allRefConcs;
                
                subIRefConc = ensemble.populations.models(modelI).metConcRef(subI) * 10^3;
            
                for subConc=substrateRange
                    subsConc(subOrderPos) = subConc ./ subIRefConc; 
                    X = [subsConc; zeros(nInhib+nActiv+nProds,1)];   

                    [v, ~, ~] = rateLawFxn(X,K);
                    vList = [vList; v];
                end

                subList = [subList; substrateRange'];
                modelList = [modelList; ones(50,1) * modelI];

            end
            write(table(modelList, subList, vList), fullfile(outputFolder, [ensemble.rxns{rxnI},'_', ensemble.mets{subI},'.csv']));
        end
    end
end

    