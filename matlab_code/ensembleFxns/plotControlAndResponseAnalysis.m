function plotControlAndResponseAnalysis(mcaResults, ensemble)
% Plots the results from the control and response analysis, four subplots 
% are generated: one for flux control coefficients, one for concentration
% control, one for enzyme response coefficients, and one for concentration
% response coefficients.
%
% To plot only some of the reactions/enzymes, the categories argument
% should be used to specify which reactions should be included in the 
% plot.
% The categories cell should look like
%
%  {reactionGroupName1, listOfReactions1; reactionGroupName2, listOfReactions2}
% 
% Example: 
%
%  {'Glycolysis',[1,20]; 'Pentose Phosphate Pathway',[25,34]}
%
%
% USAGE:
%
%    plotControlAndResponseAnalysis(mcaResults, ensemble, categories)
%
% INPUT:
%    mcaResults (struct):   MCA results
%    ensemble (struct):     model ensemble
%    categories (cell):     categories to be included in the plot
%
% .. Authors:
%       - Pedro Saa         2018 original code
%       - Marta Matos       2018 extended to plot response coefficients
%       - Marina de Leeuw   2019 added categories

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

    subplot(4,1,1)
    imagesc(mcaResults.vControlAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:numFluxes,'ytick',1:numFluxes)
    set(gca,'xticklabel',rxnNames,'yticklabel',rxnNames)
    xlabel('Reactions')
    ylabel('Reactions')
    title(['Flux control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax=subplot(4,1,2);
    colormap(ax,cmap)
    
    subplot(4,1,2)
    imagesc(mcaResults.xControlAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:numel(ix_mets),'xtick',1:numFluxes)
    set(gca,'yticklabel',metNames,'xticklabel',rxnNames)
    xlabel('Reactions')
    ylabel('Metabolites')
    title(['Concentration control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax = subplot(4,1,1);
    colormap(ax,cmap)
    colormap(cmap)
    
    
    subplot(4,1,3)
    imagesc(mcaResults.eResponseAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:size(mcaResults.eResponseAvg{ix},2),'ytick',1:numFluxes)
    set(gca,'xticklabel',mcaResults.enzNames,'yticklabel',rxnNames)
    xlabel('Enzymes')
    ylabel('Reactions')
    title(['Enzyme control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax=subplot(4,1,3);
    colormap(ax,cmap)

    subplot(4,1,4)
    imagesc(mcaResults.xResponseAvg{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:size(mcaResults.xResponseAvg{ix},2),'ytick',1:numFluxes)
    set(gca,'xticklabel',mcaResults.enzNames,'yticklabel',metNames)
    xlabel('Enzymes')
    ylabel('Metabolites')
    title(['Enzyme control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-3 3])
    ax=subplot(4,1,4);
    colormap(ax,cmap)

end


end