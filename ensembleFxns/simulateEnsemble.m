function simulationRes = simulateEnsemble(ensemble, enzymesIC, metsIC)
%
% Takes in a model ensemble, and initial conditions for enzymes and 
% metabolite concentrations and simulates all models in the ensemble.
%
%---------------- Pedro Saa UQ 2018, Marta Matos 2019 ---------------------



strucIdx = 1;
if ensemble.populations(end).strucIdx(1)==0
    ensemble.populations(end).strucIdx = ones(numel(ensemble.populations(end).strucIdx),1);
end

% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
numModels = size(ensemble.populations.models, 2);
if numModels > numel(particleIdx) 
    numModels   = numel(particleIdx);
end

% Optimization & simulation parameters
fixedExchs   = ensemble.fixedExch;
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
Sred         = ensemble.Sred;
kinInactRxns = ensemble.kinInactRxns;
subunits     = ensemble.subunits{strucIdx};

ix = 1;
simulationRes = {};

folderName =  strcat('reactions_', ensemble.description, '_', num2str(ix));
if isfile(fullfile(folderName, strcat(func2str(kineticFxn), '_ode.m')))
    odeFunction = str2func(strcat(func2str(kineticFxn), '_ode'));
else
    error(['You need a model function to be used for the model ode simulations. It should be named as ', strcat(func2str(kineticFxn), '_ode')]);
end

for jx = 1:numModels
    disp(strcat('model_', num2str(jx)));
    
    model = ensemble.populations(end).models(particleIdx(jx));

    % Simulate metabolite concentrations
    [t, y] = ode15s(@(t,y) odeFunction(y,enzymesIC,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits,0), [0,1000], metsIC);
    simulationRes{jx}.t = t;
    simulationRes{jx}.y = y;   

end

