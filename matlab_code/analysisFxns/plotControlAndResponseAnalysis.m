function plotControlAndResponseAnalysis(mcaResults, ensemble, categories, enzymeCategories)
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
%    ensemble (struct):     model ensemble, see buildEnsemble for fields description
%    categories (cell):     categories to be included in the plot
%
% .. Authors:
%       - Pedro Saa         2018 original code
%       - Marta Matos       2018 extended to plot response coefficients
%       - Marina de Leeuw   2019 added categories

% Check sampler mode to determine the numer of conditions
if ~strcmpi(ensemble.sampler,'GRASP')
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

% Checks whether any categories were defined
if (nargin < 3) || isempty(categories)
    categories = {'ControlAndResponse_reactions',[1, length(rxnNames)]};
end

if (nargin < 4) ||isempty(enzymeCategories)
    enzymeCategories = {'ControlAndResponse_enzymes',[1, length(mcaResults.enzNames)]};
end

for ix = 1:nCondition
    % Plot final results
    for j = 1:size(categories,1)
    	figure ('Name', [categories{j,1} ' , condition: ' num2str(ix)])

        subplot(2,1,1)
        imagesc(mcaResults.vControlAvg{ix}(:,categories{j,2}(1):categories{j,2}(2)))
        set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:numFluxes,'ytick',1:numFluxes)
        set(gca,'xticklabel',rxnNames(categories{j,2}(1):categories{j,2}(2)),'yticklabel',rxnNames)
        xlabel('Reactions')
        ylabel('Reactions')
        title(['Flux control coefficients condition: ',num2str(ix)])
        set(gca,'FontSize',6,'FontName','arial')
        caxis([-3 3])
        ax=subplot(2,1,2);
        colormap(ax,cmap)

        subplot(2,1,2)
        imagesc(mcaResults.xControlAvg{ix}(:,categories{j,2}(1):categories{j,2}(2)))
        set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:numel(ix_mets),'xtick',1:numFluxes)
        set(gca,'yticklabel',metNames,'xticklabel',rxnNames(categories{j,2}(1):categories{j,2}(2)))
        xlabel('Reactions')
        ylabel('Metabolites')
        title(['Concentration control coefficients condition: ',num2str(ix)])
        set(gca,'FontSize',6,'FontName','arial')
        caxis([-3 3])
        ax = subplot(2,1,1);
        colormap(ax,cmap)
        colormap(cmap)
    end
    
    for l = 1:size(enzymeCategories,1)
        figure ('Name', [enzymeCategories{l,1} ' , condition: ' num2str(ix)])
        subplot(2,1,1)
        imagesc(mcaResults.eResponseAvg{ix}(:,enzymeCategories{j,2}(1):enzymeCategories{j,2}(2)))
        set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:size(mcaResults.eResponseAvg{ix},2),'ytick',1:numFluxes)
        set(gca,'xticklabel',mcaResults.enzNames(enzymeCategories{j,2}(1):enzymeCategories{j,2}(2)),'yticklabel',rxnNames)
        xlabel('Enzymes')
        ylabel('Reactions')
        title(['Enzyme control coefficients condition: ',num2str(ix)])
        set(gca,'FontSize',6,'FontName','arial')
        caxis([-3 3])
        ax=subplot(2,1,1);
        colormap(ax,cmap)
        
        subplot(2,1,2)
        imagesc(mcaResults.xResponseAvg{ix}(:,enzymeCategories{j,2}(1):enzymeCategories{j,2}(2)))
        set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:size(mcaResults.xResponseAvg{ix},2),'ytick',1:numFluxes)
        set(gca,'xticklabel',mcaResults.enzNames(enzymeCategories{j,2}(1):enzymeCategories{j,2}(2)),'yticklabel',metNames)
        xlabel('Enzymes')
        ylabel('Metabolites')
        title(['Enzyme control coefficients condition: ',num2str(ix)])
        set(gca,'FontSize',6,'FontName','arial')
        caxis([-3 3])
        ax=subplot(2,1,2);
        colormap(ax,cmap)
    end

end


end