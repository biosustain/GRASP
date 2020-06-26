function [ensemble, models] = sampleGibbsReactionEnergies(ensemble, models, strucIdx)
% Function used to sample Gibbs energies for each reaction.
%
% The Gibbs energies for each reaction are calculated as a function of
% the standard Gibbs energies and metabolite concentrations.
% 
% Standard Gibbs energies and metabolite concentrations are found by 
% sampling their values in a thermodynamically consistent way.
% Initial boundaries for each variable and for the reaction Gibbs energy 
% are previously calculated by using TMFA.
%
% However, since the values of standard Gibbs energies and metabolites 
% are not independent of each other, we sample one variable at a time
% within the  boundaries previously calculated. However, after sampling 
% each variable, the boundaries for the remaining variables change, so  
% before sampling a variable its boundaries are updated by solving a linear 
% program to find the minimum and maximum possible values of that variable.
%   
% The constraints for the linear program are the following:
%
% .. math::
%
%       \Delta G_{min} \leq \Delta G <= \Delta G_{max}
%                      
% which is equivalent to:
%
% .. math::
%
%       \Delta G_{min} <= \Delta G^0 + RT S' ln(x) <= \Delta G_{max}
% 
% where S is the stoichiometric matrix, R is the gas constant, T is 
% temperature, x are metabolite concentrations, and G^0 are standard Gibbs
% energies. The metabolite concentrations and standard Gibbs energies are
% the variables to be sampled and are bounded:
%
% .. math::
%
%       \Delta G^0_{min} <= \Delta G^0 <= \Delta G^0_{max} 
%
%       x_{min} <= x <= x_{max}
%
% After calculating the the new bounds for a given variable, its lower and
% upper bounds are set to the sampled value and the linear program is
% updated. This procedure is repeated for each variable.
%
%
% Both Gurobi and linprog can be used to solv ethe linear programs.
%
%
% USAGE:
%
%    [ensemble, models] = sampleGeneralReversibilities(ensemble, models, strucIdx)
%
% INPUT:
%    ensemble (struct):	  model ensemble, see buildEnsemble for fields description
%    models (struct):     model, see initialSampler for fields description
%    strucIdx (int):      number of the model structure considered
%
% OUTPUT:
%    ensemble (struct):	  model ensemble, see buildEnsemble for fields description
%    models (struct):     model data structure, see initialSampler for fields description
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
    
    A = sparse([-eye(nRxns), -RT*ensemble.Sthermo';...
                eye(nRxns), RT*ensemble.Sthermo']);

    b = [-ensemble.gibbsRanges((ensemble.thermoActive),1);... 
         ensemble.gibbsRanges((ensemble.thermoActive),2)];

    f = zeros(nRxns+nMets,1);
    
    lb = [ensemble.G0Ranges(ensemble.thermoActive,1); log(ensemble.metRanges(:,1))];
    ub = [ensemble.G0Ranges(ensemble.thermoActive,2); log(ensemble.metRanges(:,2))];
       
    
    % Set up LP model for gurobi if this is the chosen solver.
    if strcmp(ensemble.LPSolver, 'gurobi')        
        gurobiModel.A = A;

        sense = zeros(1,nRxns*2);
        sense(1:nRxns) = '<';
        sense((nRxns + 1):nRxns*2) = '<';
        gurobiModel.sense = char(sense);

        gurobiModel.rhs = b;
                       
        gurobiModel.lb = lb;
        gurobiModel.ub = ub;

        gurobiModel.vtype = 'C';
        gurobiModel.obj = f; 
    end

    varList = 1:nMets+nRxns;

    % Sample standard Gibbs energies
    offset = 0;
    rxnList = randperm(nRxns);
    
    if strcmp(ensemble.LPSolver, 'gurobi')
        [lb, ub, varList] = sampleVariablesGurobi(gurobiModel, lb, ub, varList, rxnList, G0Factor, offset);
    elseif strcmp(ensemble.LPSolver, 'linprog')
        [lb, ub, varList] = sampleVariablesLinprog(f, A, b, lb, ub, varList, rxnList, G0Factor, offset);
    end    

    % Sample measured metabolites
    if numel(ensemble.measuredMets) > 0
        offset = nRxns;
        metList = ensemble.measuredMets(randperm(numel(ensemble.measuredMets)));   

        if strcmp(ensemble.LPSolver, 'gurobi')
            [lb, ub, varList] = sampleVariablesGurobi(gurobiModel, lb, ub, varList, metList, metsFactor, offset);
        elseif strcmp(ensemble.LPSolver, 'linprog')
            [lb, ub, varList] = sampleVariablesLinprog(f, A, b, lb, ub, varList, metList, metsFactor, offset);
        end

    end

    % Sample not measured metabolites
    if numel(ensemble.measuredMets) <= nMets
        offset = nRxns;
        allMets = 1:nMets;
        notMeasMets = allMets(~ismember(allMets, ensemble.measuredMets));
        metList = notMeasMets(randperm(numel(notMeasMets)));    

        if strcmp(ensemble.LPSolver, 'gurobi')
            [lb, ub, varList] = sampleVariablesGurobi(gurobiModel, lb, ub, varList, metList, metsFactor, offset);
        elseif strcmp(ensemble.LPSolver, 'linprog')
            [lb, ub, varList] = sampleVariablesLinprog(f, A, b, lb, ub, varList, metList, metsFactor, offset);
        end

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


function [lb, ub, varList] = sampleVariablesGurobi(gurobiModel, lb, ub, varList, varSubList, varFactor, offset)
%
% Get new bounds for each variable to be sampled, by minimizing and
% maximizing its value and then sample a value within those bounds.
% At the end update the LP's lower and upper bounds.%
%

gurobiModel.lb = lb;
gurobiModel.ub = ub;
params.outputflag = 0;
params.OptimalityTol = 1e-6;
params.FeasibilityTol = 1e-6;

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
    
    if strcmp(solmin.status, 'INFEASIBLE') || strcmp(solmax.status, 'INFEASIBLE')
        error('The linear program is infeasible. Something went wrong with sampling Gibbs energies.');
    end
    
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


function [lb, ub, varList] = sampleVariablesLinprog(f, A, b, lb, ub, varList, varSubList, varFactor, offset)
%
% Get new bounds for each variable to be sampled, by minimizing and
% maximizing its value and then sample a value within those bounds.
% At the end update the LP's lower and upper bounds.%
%

nVars = numel(varSubList);
counter = 1;

options =  optimoptions(@linprog, 'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, 'Display', 'off');

while max(abs(ub - lb)) > 10^-5 && counter <= nVars
    
    varI = varSubList(counter);
    
    % Find new bounds for the variable
    f(offset+varI)    = 1;
    [x, minfval]      = linprog(f, A, b, [], [], lb, ub, options);
    f(offset+varI)    = -1;
    [x, maxfval]      = linprog(f, A, b, [], [], lb, ub, options);
    maxfval = -maxfval; 
    
    if isempty(minfval) || isempty(maxfval)
        error('The linear program is infeasible. Something went wrong with sampling Gibbs energies.');
    end    
    
    % Sample value within the new bounds
    varValue = varFactor(varI).*maxfval + (1-varFactor(varI)).*minfval;

    % Update variable bounds
    lb(offset+varI) = varValue;
    ub(offset+varI) = varValue;

    % Re-set objective function
    f(offset+varI)= 0;            

    % Remove variable from list of variables to be sampled.
    varList(find(varList==(offset+varI))) = [];
    
    counter = counter + 1;

end
end