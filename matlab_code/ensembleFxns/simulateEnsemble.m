function simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, interruptTime)
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

currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');  
folderName =  fullfile(currentPath{1}, '..', '..', 'reactions', strcat(ensemble.description, '_', num2str(ix)));
if isfile(fullfile(folderName, strcat(func2str(kineticFxn), '_ode.m')))
    odeFunction = str2func(strcat(func2str(kineticFxn), '_ode'));
else
    error(['You need a model function to be used for the model ode simulations. It should be named as ', strcat(func2str(kineticFxn), '_ode')]);
end

simulationRes = cell(1, numModels);

disp ('Simulating models.');

for jx = 1:numModels
   
    model = ensemble.populations(end).models(particleIdx(jx));
    metConcRef = model.metConcRef(ensemble.metsBalanced);
    
    outputFun= @(t,y,flag)interuptFun(t,y,flag,interruptTime);
    opts = odeset('RelTol',1e-13,'OutputFcn',outputFun);

    try
        % Simulate metabolite concentrations
        [t, y] = ode15s(@(t,y) odeFunction(y,enzymesIC,metConcRef,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits), [0,finalTime], metsIC, opts);

        simulationRes{jx}.t = t;
        simulationRes{jx}.conc = y;   
        simulationRes{jx}.flux = calculateFluxes(t,y,enzymesIC,kineticFxn,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits);   
        simulationRes{jx}.flux = simulationRes{jx}.flux ./ ensemble.fluxRef';
        
    catch ME
        if strcmp(ME.identifier,'interuptFun:Interupt')
            disp(ME.message);
        else
            rethrow(ME); % It's possible the error was due to something else
        end
    end

end

end


function status = interuptFun(t,y,flag,interruptTime)   
%
% Interrupts ODE solver if it takes more than interruptTime (in seconds);
%

persistent INIT_TIME;
status = 0;

switch(flag)
    case 'init'
        INIT_TIME = tic;
    case 'done'
        clear INIT_TIME;
    otherwise
        elapsedTime = toc(INIT_TIME);
        if elapsedTime > interruptTime
            clear INIT_TIME;
            error('interuptFun:Interupt',...
                 ['Interupted integration. Elapsed time is ' sprintf('%.6f',elapsedTime) ' seconds.']);
        end

end
end

function flux = calculateFluxes(timePoints,metConcs,enzymesIC,kineticFxn,model,fixedExchs,Sred,kinInactRxns,subunits)
    
flux = zeros(numel(timePoints), size(Sred,2));

for t=1:numel(timePoints)
    x = [metConcs(t,:)'; enzymesIC];
    flux(t,:) = feval(kineticFxn,x,model,fixedExchs,Sred,kinInactRxns,subunits,0);
end
end
