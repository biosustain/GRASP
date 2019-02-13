function plotControlAnalysis(mcaResults, ensemble)
%---------------- Pedro Saa UQ 2018----------------------------------------


% Check sampler mode to determine the numer of conditions
if ~strcmpi(ensemble.sampler,'ORACLE')
    nCondition   = size(ensemble.expFluxes,2)+1;
else
    nCondition = 1;
end


% Define colormap for the heatmap
try
    load('cmap_rgb.mat')
catch
    disp('No predifined colormap')
    cmap = colormap(jet);
end

% Optimization & simulation parameters
freeVars     = numel(ensemble.freeVars);
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsActive);
metNames     = ensemble.mets(ensemble.metsActive);
rxnNames     = ensemble.rxns;


for ix = 1:nCondition
    % Plot final results
    figure (ix)
    subplot(2,1,1)
    imagesc(mcaResults.xControlAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:numel(ix_mets),'xtick',1:numFluxes)
    set(gca,'yticklabel',metNames,'xticklabel',rxnNames)
    ylabel('Metabolites')
    xlabel('Reactions')
    title(['Concentration control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax = subplot(2,1,1);
    colormap(ax,cmap)
    colormap(cmap)
    
    subplot(2,1,2)
    imagesc(mcaResults.vControlAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:numFluxes,'ytick',1:numFluxes)
    set(gca,'xticklabel',rxnNames,'yticklabel',rxnNames)
    xlabel('Reactions')
    ylabel('Reactions')
    title(['Flux control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax=subplot(2,1,2);
    colormap(ax,cmap)
    

end


end