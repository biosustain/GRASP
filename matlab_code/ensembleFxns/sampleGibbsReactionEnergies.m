function [ensemble, models] = sampleGibbsReactionEnergies(ensemble, models, strucIdx)
% Function used to sample Gibbs energies for each reaction.
%
% A = G0 + RT S' ln(x)
%
% variables: G0, X=ln(x)
%
% The method is going to change, so no need to write a good description
%
%
% USAGE:
%
%    [ensemble, models] = sampleGeneralReversibilities(ensemble, models, strucIdx)
%
% INPUT:
%    ensemble (struct):	  model ensemble
%    models (struct):     model
%    strucIdx (int):      number of the model structure considered
%
% OUTPUT:
%    ensemble (struct):	  model ensemble
%    models (struct):     model data structure 
%
% .. Authors:
%       - Marta Matos       2019 original code

RT = 8.314*298.15/1e3; 

if isfield(ensemble,'G0Ranges') == 1
    nMets = size(ensemble.metRanges, 1);
    nRxns = size(ensemble.G0Ranges(ensemble.thermoActive), 1);


    % Sample metabolite concentrations
    metsFactor = mvnrnd(zeros(nMets, 1), eye(nMets))';
    metsFactor = exp(metsFactor)./(1 + exp(metsFactor));

    % Sample standard Gibbs energies
    G0Factor = mvnrnd(zeros(nRxns, 1),eye(nRxns))';
    G0Factor = exp(G0Factor)./(1 + exp(G0Factor));


    % Set up gurobi model
    gurobiModel.A = sparse([eye(nRxns), RT*ensemble.Sthermo';...
                            eye(nRxns), RT*ensemble.Sthermo']);

    sense = zeros(1,nRxns*2);
    sense(1:nRxns) = '>';
    sense((nRxns + 1):nRxns*2) = '<';
    gurobiModel.sense = char(sense);


    gurobiModel.rhs = [ensemble.gibbsRanges((ensemble.thermoActive),1);... 
                       ensemble.gibbsRanges((ensemble.thermoActive),2)];


    lb = [ensemble.G0Ranges(ensemble.thermoActive,1); log(ensemble.metRanges(:,1))];
    ub = [ensemble.G0Ranges(ensemble.thermoActive,2); log(ensemble.metRanges(:,2))];
    gurobiModel.lb = lb;
    gurobiModel.ub = ub;

    gurobiModel.vtype = 'C';
    gurobiModel.obj = zeros(nRxns+nMets,1);

    varList = 1:nMets+nRxns;

    % Sample standard Gibbs energies
    offset = 0;
    rxnList = randperm(nRxns);
    [lb, ub, varList] = sampleVariables(gurobiModel, lb, ub, varList, rxnList, G0Factor, offset);

    % Sample measured metabolites
    if numel(ensemble.measuredMets) > 0
        offset = nRxns;
        metList = ensemble.measuredMets(randperm(numel(ensemble.measuredMets)));   
        [lb, ub, varList] = sampleVariables(gurobiModel, lb, ub, varList, metList, metsFactor, offset);
    end

    % Sample not measured metabolites
    if numel(ensemble.measuredMets) <= nMets
        offset = nRxns;
        allMets = 1:nMets;
        notMeasMets = allMets(~ismember(allMets, ensemble.measuredMets));
        metList = notMeasMets(randperm(numel(notMeasMets)));    
        [lb, ub, varList] = sampleVariables(gurobiModel, lb, ub, varList, metList, metsFactor, offset);
    end

    assert(max(abs(ub-lb)) < 10^-5, 'dG values are not fully determined');
    
    % Calculate Gibbs energies and save sampled metabolite reference concentrations.
    ensemble.gibbsTemp = -1e2*ones(size(ensemble.G0Ranges, 1), 1);
    ensemble.gibbsTemp(ensemble.thermoActive) = lb(1:nRxns) + RT*ensemble.Sthermo'*lb((nRxns+1):end);
    models(1).metConcRef = exp(lb((nRxns+1):end));
    models(1).gibbsTemp = ensemble.gibbsTemp;


else
    disp([newline,'Warning: sampling Gibbs energies directly might result in thermodynamically infeasible models!', newline]);
    
    gibbsFactor = mvnrnd(ensemble.populations(1).probParams(strucIdx).muGibbsFactor,ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor)';
    gibbsFactor = exp(gibbsFactor)./(1 + exp(gibbsFactor));
    gibbsEnergy = gibbsFactor.*ensemble.gibbsRanges(ensemble.thermoActive,2) + (1-gibbsFactor).*ensemble.gibbsRanges(ensemble.thermoActive,1);
    models(1).gibbsFactor = gibbsFactor;
    
    ensemble.gibbsTemp = zeros(size(ensemble.thermoActive));
    thermoCounter = 1;
    for activRxnIdx = 1:numel(ensemble.kinActRxns)        

        if ~(strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
             strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange'))

            % Determine gibbs free energy of reaction
            if ismember(ensemble.kinActRxns(activRxnIdx),ensemble.thermoActive)
                ensemble.gibbsTemp(ensemble.kinActRxns(activRxnIdx))     = gibbsEnergy(thermoCounter);
                thermoCounter = thermoCounter+1;
            elseif (ensemble.gibbsRanges(activRxnIdx,1)==ensemble.gibbsRanges(activRxnIdx,2))
                ensemble.gibbsTemp(ensemble.kinActRxns(activRxnIdx))     = ensemble.gibbsRanges(activRxnIdx,1);
            end

        end
    end
end

end


function [lb, ub, varList] = sampleVariables(gurobiModel, lb, ub, varList, varSubList, varFactor, offset)
%
% Get new bounds for each variable to be sampled, by minimizing and
% maximizing its value and then sample a value within those bounds.
% At the end update the LP's lower and upper bounds.%
%

gurobiModel.lb = lb;
gurobiModel.ub = ub;
params.outputflag = 0;

nVars = numel(varSubList);
counter = 1;

while max(abs(ub - lb)) > 10^-5 && counter <= nVars
    
    varI = varSubList(counter);

    % Find new bounds for the variable
    gurobiModel.obj(offset+varI)    = 1;
    gurobiModel.modelsense = 'min';
    solmin           = gurobi(gurobiModel,params);
    gurobiModel.modelsense = 'max';
    solmax           = gurobi(gurobiModel,params);

    % Sample value within the new bounds
    varValue = varFactor(varI).*solmax.objval + (1-varFactor(varI)).*solmin.objval;

    % Update variable bounds
    lb(offset+varI) = varValue;
    ub(offset+varI) = varValue;
    gurobiModel.lb = lb;
    gurobiModel.ub = ub;

    % Re-set objective function
    gurobiModel.obj(offset+varI)= 0;            

    % Remove variable from list of variables to be sampled.
    varList(find(varList==(offset+varI))) = [];
    
    counter = counter + 1;

end

end