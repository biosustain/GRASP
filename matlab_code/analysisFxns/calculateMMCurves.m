function calculateMMCurves(outputFolder, ensemble, numModels, strucIdx)

% in mmol/L
%substrateRange = logspace(-6, 4);
%saturatingConc = 10^4;

saturatingConc = 10^2*10^6;
substrateRange = logspace(-15,6);

%for rxnI=1:numel(ensemble.rxns)
for rxnI=32:32
    if sum(~ismember({'massAction', 'fixedExchange', 'freeExchange', 'diffusion'}, ensemble.rxnMechanisms{strucIdx}{rxnI})) == 4 
       
        rateLawFxn = str2func([ensemble.rxns{rxnI}, num2str(strucIdx)]);

        stoicSubsInd = find(ensemble.S(:, rxnI) < 0);
        stoicProdInd = find(ensemble.S(:, rxnI) > 0);
        
        nSubs = numel(ensemble.subOrder{1}{rxnI});

        for i=1:numel(stoicSubsInd)
            subI = stoicSubsInd(i);
            
            subList = [];
            vList = [];
            modelList = [];
            
            % catch: concentrations of substrates for promiscuous
            % reactions must be zero instead of saturating.
            subsConc = zeros(nSubs, 1);
            coSubsInd = find(ismember(ensemble.subOrder{1}{rxnI}, ensemble.mets(stoicSubsInd)));
            
            % met position in ensemble.subOrder
            subOrderPos = find(ismember(ensemble.subOrder{1}{rxnI}, ensemble.mets(subI)));
            
            % ensemble.Order met position in ensemble.mets
            subOrderInd = [];
            for entry=1:numel(ensemble.subOrder{1}{rxnI})
                subOrderInd = [subOrderInd, find(ismember(ensemble.mets, ensemble.subOrder{1}{rxnI}{entry}))];
            end
            
            for modelI=1:2    

                K = ensemble.populations.models(modelI).rxnParams(rxnI).kineticParams;
                allRefConcs = ensemble.populations.models(1).metConcRef(subOrderInd) * 10^6;
                                
                subsConc(coSubsInd) = saturatingConc ./ allRefConcs;
                
                subIRefConc = ensemble.populations.models(1).metConcRef(subI) * 10^6;
            
                for subConc=substrateRange
                    subsConc(subOrderPos) = subConc ./ subIRefConc; 
                    X = [subsConc; zeros(numel(stoicProdInd),1)];   

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

    