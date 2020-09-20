function simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime, numModels, numCores)
% Takes in a model ensemble and initial conditions for enzyme and 
% metabolite concentrations and simulates all models in the ensemble.
%
% For metabolites, the initial conditions can be given either in terms of
% relative or absolute concentrations (mol/L) by setting *metsAbsOrRel* to
% 'rel' or 'abs', respectively.
%
% Metabolite concentrations are always returned as relative values while
% reaction fluxes are always returned as absolute values.
%
% If the simulation of a given model takes longer than the specified 
% *interrupTime*, then it is interrupted and no simulation results are  
% saved for that model.
%
% Note that for each model a different set of time points and respective
% concentrations/fluxes will be returned, since the solver's step size is
% adaptive, i.e., not constant.
%
%
% USAGE:
%
%    simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, interruptTime)
%
% INPUT:
%    ensemble (struct):           model ensemble, see buildEnsemble for fields description
%    finalTime (double):          simulation time
%    enzymesIC (double vector):	  initial conditions for enzyme concentrations
%    metsIC (double vector):      initial conditions for metabolite concentrations
%    metsAbsOrRel (char):         specify whether metabolite initial conditions are relative or absolute concentrations. Accepted values for this variable are 'abs' and 'rel'
%    interruptTime (double)       maximum time for each simulation, given in seconds
%    numModels (int):             how many models should be simulated. This number should be lower than the number of models in the ensemble.
%
% OPTIONAL INPUT:
%    numCores (int):              number of cores to be used to run the code, default is 1.
%
% OUTPUT:
%    simulationRes (struct):  simulation results
%
%               * t (*cell*)      : time points in each model simulation
%               * conc (*cell*)   : relative concentrations for each time point and model simulation
%               * flux (*cell*)   : absolute fluxes for each time point and model simulation
%
% .. Authors:
%       - Marta Matos       2019 original code

if ~strcmp(metsAbsOrRel, 'rel') && ~strcmp(metsAbsOrRel, 'abs')
    error('The value of the variable metsAbsOrRel must be either "rel" or "abs".');
end

if nargin < 8
    numCores = 1;
end

strucIdx = 1;
if ensemble.populations(end).strucIdx(1)==0
    ensemble.populations(end).strucIdx = ones(numel(ensemble.populations(end).strucIdx),1);
end

% Add kinetic fxns to the path
addKineticFxnsToPath(ensemble);

% Find particles of the appropriate structure
particleIdx = find(ensemble.populations(end).strucIdx==strucIdx);
if numModels > numel(particleIdx) 
    numModels   = numel(particleIdx);
end

% Optimization & simulation parameters
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
Sred         = ensemble.Sred;
kinInactRxns = ensemble.kinInactRxns;
subunits     = ensemble.subunits{strucIdx};
xconstIC     = ones(numel(ensemble.metsFixed), 1);
xvarIC       = ones(numel(ensemble.metsActive), 1);

nfreeVars = numel(ensemble.freeVars);
enzIC = ones(nfreeVars-numel(ensemble.metsActive) ,1);

if numel(metsIC) > 0
    [xvarIC, xconstIC] = initializeMetIC(ensemble, xvarIC, xconstIC, metsIC);
end

if numel(enzymesIC) > 0
    enzIC = initializeEnzIC(ensemble, enzIC, enzymesIC);
end

ix = 1;

currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');  
folderName =  fullfile(currentPath{1}, '..', '..', 'reactions', strcat(ensemble.description, '_', num2str(ix)));
if isfile(fullfile(folderName, strcat(func2str(kineticFxn), '_ode.m')))
    odeFunction = str2func(strcat(func2str(kineticFxn), '_ode'));
else
    error(['You need a model function to be used for the model ode simulations. It should be named as ', strcat(func2str(kineticFxn), '_ode')]);
end

timePoints = logspace(-10, log10(finalTime), 100);

xvarICabs = xvarIC;  % dirty trick to make parfor work -.-
xconstICabs = xconstIC;  % dirty trick to make parfor work -.-
simulationRes = cell(1, numModels);

disp ('Simulating models.');

parpool(numCores);																								% Initiate parallel pool and run parallel foor loop
parfor jx = 1:numModels
    disp(['Model: ', num2str(jx)]);
    
    fixedExchs = ensemble.populations(end).models(jx).fixedExch;
    
    model = ensemble.populations(end).models(particleIdx(jx));
    metActiveConcRef = model.metConcRef(ensemble.metsActive);
    metFixedConcRef = model.metConcRef(ensemble.metsFixed);
    
    if strcmp(metsAbsOrRel, 'abs')
        perturbInd = find(xvarICabs ~= 1);
        xvarICtemp = ones(size(xvarICabs));
        xvarICtemp(perturbInd) = xvarICabs(perturbInd) ./ metActiveConcRef(perturbInd); 
    
        perturbInd = find(xconstICabs ~= 1);
        xconstICtemp = ones(size(xconstICabs));
        xconstICtemp(perturbInd) = xconstICabs(perturbInd) ./ metFixedConcRef(perturbInd);
        
    elseif strcmp(metsAbsOrRel, 'rel')
        xvarICtemp = xvarIC;   
        xconstICtemp = xconstIC;
    end
    
    outputFun= @(t,y,flag)interuptFun(t,y,flag,interruptTime);
    opts = odeset('RelTol',1e-13,'OutputFcn',outputFun);

    try
        % Simulate metabolite concentrations
        [t, y] = ode15s(@(t,y) odeFunction(y,enzIC,metActiveConcRef,xconstICtemp,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits), timePoints, xvarICtemp, opts);

        simulationRes{jx}.t = t;
        simulationRes{jx}.conc = y;   
        simulationRes{jx}.flux = calculateFluxes(t,y,enzIC,kineticFxn,xconstICtemp,model,fixedExchs(:,ix),Sred,kinInactRxns,subunits);   
        
    catch ME
        if strcmp(ME.identifier,'interuptFun:Interupt')
            disp(ME.message);
        else
            rethrow(ME); % It's possible the error was due to something else
        end
    end

end
delete(gcp);

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


function flux = calculateFluxes(timePoints,metConcs,enzymesIC,kineticFxn,xconst,model,fixedExchs,Sred,kinInactRxns,subunits)
    
flux = zeros(numel(timePoints), size(Sred,2));

for t=1:numel(timePoints)
    x = [metConcs(t,:)'; enzymesIC];
    flux(t,:) = feval(kineticFxn,x,xconst,model,fixedExchs,Sred,kinInactRxns,subunits,0);
end

end



function [xvar, xconst] = initializeMetIC(ensemble, xvar, xconst, metsIC)
%
% Initializes xvar and xconst with the metabolite concentration values 
% given by the user.
%
   
for metI=1:numel(metsIC)
    metInd = find(ismember(ensemble.mets, ['m_', metsIC{metI}{1}]) == 1);
    
    if numel(metInd) == 0
        error(strcat("The metabolite ", metsIC{metI}{1}, " doesn't exist in the model"));
    end
    
    activeMetInd = find(ensemble.metsActive == metInd);
    
    if numel(activeMetInd) == 0
        fixedMetInd = find(ensemble.metsFixed == metInd);
        xconst(fixedMetInd) = metsIC{metI}{2};
    else
        xvar(activeMetInd) = metsIC{metI}{2};
    end
end

end


function enzIC = initializeEnzIC(ensemble, enzIC, enzymesIC)
%
% Initializes enzIC with the enzyme concentrations given by the user.
%
    
for enzI=1:numel(enzymesIC)
    enzInd = find(ismember(ensemble.rxns, ['r_', enzymesIC{enzI}{1}]) == 1);
    
    if numel(enzInd) == 0
        error(strcat("The reaction ", enzymesIC{enzI}{1}, " doesn't exist in the model"));
    end
    
    enzIC(enzInd) = enzymesIC{enzI}{2};
end

end
