function controlAnalysis(ensemble,strucIdx)
%--------------------------- Pedro Saa UQ 2018 ----------------------------
% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
numModels   = numel(particleIdx);

% Optimization & simulation parameters
nCondition   = size(ensemble.expFluxes,2)+1;
fixedExchs   = ensemble.fixedExch;
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = size(ensemble.populations(end).xopt{1},1);
Sred         = ensemble.Sred;
kinInactRxns = ensemble.kinInactRxns;
subunits     = ensemble.subunits{strucIdx};
numFluxes    = size(ensemble.populations(end).simFluxes{1},1);
ix_enz       = ix_mets(end)+1:freeVars;
metNames     = ensemble.mets(ensemble.metsActive);
rxnNames     = ensemble.rxns;

% Define colormap for the heatmap
try
    load('cmap_rgb.mat')
catch
    disp('No predifined colormap')
    cmap = colormap(jet);
end

% Main loop
hstep = 1e-10;              % Step size for control coefficient computations
for ix = 1:nCondition
    xControl{ix} = 0;
    eControl{ix} = 0;
    xcounter{ix} = 0;
    ecounter{ix} = 0;
    for jx = 1:numModels
        model = ensemble.populations(end).models(particleIdx(jx));
        if ix == 1
            xopt = ones(freeVars,1);
            vref = feval(kineticFxn,xopt,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0);
        else
            xopt = ensemble.populations(end).xopt{particleIdx(jx)}(:,ix-1);
            vref = ensemble.populations(end).simFluxes{particleIdx(jx)}(:,ix-1);
        end
        
        % Define reference state
        xref = xopt(ix_mets);
        Eref = xopt(ix_enz);

        % Define step length
        hstep_x = hstep*xref;
        xmets   = repmat(xref,1,numel(xref)) + 1i*diag(hstep_x);
        xenz    = repmat(Eref,1,numel(xref));
        xstep   = [xmets;xenz];
        
        % Simulate flux
        simFlux = feval(kineticFxn,xstep,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0);
        
        % Compute elasticiy matrix
        E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))';
        
        % Compute control coefficients
        C_x_abs   = -(pinv(Sred*E_x_abs))*Sred;
        C_x       = diag(xref.^-1)*C_x_abs*diag(vref);
        E_x_nor   = diag(vref.^-1)*E_x_abs*diag(xref);
        C_e       = eye(numel(vref)) + E_x_nor*C_x;
        
        % Save control coefficients only if the result is accurate
        if all(abs(sum(C_x,2))<1e-5)
            xControl{ix} = xControl{ix} + C_x;
            xcounter{ix} = xcounter{ix} + 1;
        end
        if all(abs(sum(C_e,2))-1<1e-5)
            eControl{ix} = eControl{ix} + C_e;
            ecounter{ix} = ecounter{ix} + 1;
        end
    end
    
    % Determine expectancy for control coefficients
    xControl{ix} = xControl{ix}/xcounter{ix};
    eControl{ix} = eControl{ix}/ecounter{ix};
    
    % Plot final results
    figure (ix)
    subplot(2,1,1)
    imagesc(xControl{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'ytick',1:numel(ix_mets),'xtick',1:numFluxes)
    set(gca,'yticklabel',metNames,'xticklabel',rxnNames)    
    ylabel('Metabolites')
    xlabel('Reactions')
    title(['Concentration control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-2.5 2.5])
    ax = subplot(2,1,1);
    colormap(ax,cmap)
    colormap(cmap)
    subplot(2,1,2)
    imagesc(eControl{ix})
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',1:numFluxes,'ytick',1:numFluxes)
    set(gca,'xticklabel',rxnNames,'yticklabel',rxnNames)
    xlabel('Enzymes')
    ylabel('Reactions')
    title(['Flux control coefficients condition: ',num2str(ix)])
    set(gca,'FontSize',6,'FontName','arial')
    caxis([-2.5 2.5])
    ax=subplot(2,1,2);
    colormap(ax,cmap)
end

% Save final results
save('control_results','xControl','eControl','xcounter','ecounter')
