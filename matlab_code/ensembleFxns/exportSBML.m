function exportSBML(ensemble, modelI)

    SBmodel = sbiomodel(ensemble.description);
    SBmodel = addRefConc(SBmodel, ensemble, modelI);
    fluxVector = cell(numel(ensemble.rxns), 1);

    currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
    reactionsFolder = fullfile(currentPath{1}, '..', '..',  'reactions', ...
                             strcat(ensemble.description, '_1'));

    for rxnI=1:numel(ensemble.rxns)
        Kpos = {};
        Kneg = {};
        clearvars fluxEq freeEnz
        disp(rxnI);
               
        rxnName = ensemble.rxns{rxnI};
        catalyticFxn = fullfile(reactionsFolder, [rxnName, '1Catalytic.m']);
        fullFxn = fullfile(reactionsFolder, [rxnName, '1.m']);   
        
        if isfile(catalyticFxn)
            
            if strcmp(ensemble.rxnMechanisms{1}{rxnI}, 'massAction')
                fluxEq = parseMassActionFunction(ensemble, rxnI);  
                freeEnz = 1;
            else
                [fluxEq, freeEnz] = parseCatalyticFunction(catalyticFxn, ensemble, rxnI);  
            end
            
            [fluxEq, Kpos, Kneg] = parseAllostericFunction(ensemble, rxnI, fluxEq, freeEnz);   
            
        else
            
            if strcmp(ensemble.rxnMechanisms{1}{rxnI}, 'massAction')                
                fluxEq = parseMassActionFunction(ensemble, rxnI);
            else
                [fluxEq, freeEnz] = parseCatalyticFunction(fullFxn, ensemble, rxnI); 
            end
            
        end
        
        [SBmodel, fluxEqSub] = addReactionToSBModel(SBmodel, ensemble, fluxEq, modelI, rxnI, Kpos, Kneg, catalyticFxn);        
        fluxVector{rxnI} = str2sym(fluxEqSub);

    end
    
    SBmodel = setSpeciesInitialConcentrations(SBmodel, ensemble);
    SBmodel = setConstSpecies(SBmodel, ensemble);
    SBmodel = addODEsToSBModel(SBmodel, ensemble, fluxVector);

    sbmlexport(SBmodel,'/home/mrama/Desktop/test_1.xml');
    
    postProcessSBMLFile('/home/mrama/Desktop/test_1.xml', ensemble);

end


function fluxEq = parseMassActionFunction(ensemble, rxnI)
    syms k01 k02
    
    substrates = ensemble.mets(find(ensemble.S(:, rxnI) < 0));
    stoicCoefSub = ensemble.S(find(ensemble.S(:, rxnI) < 0), rxnI);
    forRate = '';
    for subI=1:numel(substrates)
        forRate = [forRate, ' * ', substrates{subI}, '^', num2str(abs(stoicCoefSub(subI)))];
    end
    
    
    products = ensemble.mets(find(ensemble.S(:, rxnI) > 0));
    stoicCoefProd = ensemble.S(find(ensemble.S(:, rxnI) > 0), rxnI);
    revRate = '';
    for prodI=1:numel(products)
        revRate = [revRate, ' * ', products{prodI}, '^', num2str(stoicCoefProd(prodI))];
    end

    fluxEq = k01 * str2sym(forRate(4:end)) - k02 * str2sym(revRate(4:end));

end


function [fluxEq, freeEnz] = parseCatalyticFunction(fullFxn, ensemble, rxnI)

    varList = '';
    functionContents = '';

    varList = getRateConstNames(ensemble, rxnI, varList);

    if isfile(fullFxn)
        fileIn = fopen(fullFxn);
        line = fgets(fileIn);
        
        while line ~= -1
            [functionContents, varList] = parseLine(line, ensemble, rxnI, functionContents, varList); 
            line = fgets(fileIn);
        end
        
        functionContents = [functionContents, '\nv = simplify(v);\n'];
        functionContents = [functionContents, 'E1 = simplify(E1);\n'];
        
        fclose(fileIn);
    end

    currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');

    writeToFile(currentPath, functionContents, varList);

    [fluxEq, freeEnz]  = runFluxEq(currentPath);
    
end


function varList = getRateConstNames(ensemble, rxnI, varList)

    nRates = numel(ensemble.populations.models(1).rxnParams(rxnI).kineticParams);
    rateI = 1;

    while rateI < 10 && rateI <= nRates
        varList = [varList, [' k0', num2str(rateI)]];
        rateI = rateI + 1;
    end
    
    while rateI <= nRates
        varList = [varList, [' k', num2str(rateI)]];
        rateI = rateI + 1;
    end

end


function [functionContents, varList] = parseLine(line, ensemble, rxnI, functionContents, varList)

    Eline = numel(find(strfind(line, 'E') == 1)) == 1 && numel(strfind(line, 'E')) >= 1;
    Dline = all(strfind(line, 'Denominator') == 1) && numel(strfind(line, 'Denominator')) >= 1;
    vline = all(strfind(line, 'v') == 1) && numel(strfind(line, 'v')) >= 1;

    if strfind(line, 'A =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'A', ensemble.subOrder{1,1}{rxnI}{1});

    elseif strfind(line, 'B =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'B', ensemble.subOrder{1,1}{rxnI}{2});

    elseif strfind(line, 'C =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'C', ensemble.subOrder{1,1}{rxnI}{3});

    elseif strfind(line, 'D =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'D', ensemble.subOrder{1,1}{rxnI}{4});

    elseif strfind(line, 'E =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'E', ensemble.subOrder{1,1}{rxnI}{5});

    elseif strfind(line, 'F = ') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'F', ensemble.subOrder{1,1}{rxnI}{6});

    elseif strfind(line, 'P =') == 1    
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'P', ensemble.prodOrder{1,1}{rxnI}{1});

    elseif strfind(line, 'Q =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'Q', ensemble.prodOrder{1,1}{rxnI}{2});

    elseif strfind(line, 'R =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'R', ensemble.prodOrder{1,1}{rxnI}{3});

    elseif strfind(line, 'S =') == 1
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                'S', ensemble.prodOrder{1,1}{rxnI}{4});

    elseif strfind(line, 'I') == 1
        if strcmp(line(2), ' ') == 1
            inhibNum = 1;
        else
            inhibNum = str2num(line(2));
        end

        eqPos = strfind(line, '=');
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                line(1:eqPos-2), ensemble.inhibitors{1,1}{rxnI}{inhibNum});

    elseif strfind(line, 'Z') == 1
        if strcmp(line(2), ' ') == 1
            actNum = 1;
        else
            actNum = str2num(line(2));
        end

        eqPos = strfind(line, '=');
        [functionContents, varList] = parseMets(functionContents, varList, ...
                                                line(1:eqPos-2), ensemble.activators{1,1}{rxnI}{actNum});

    elseif Eline || Dline || vline
        eqPos = strfind(line, '=');
        functionContents = [functionContents, line];
        varList = [varList, ' ', line(1:eqPos-2)];                                    
    end  

end


function [functionContents, varList] = parseMets(functionContents, varList, lhs, rhs)

    functionContents = [functionContents, [lhs, ' = ', rhs, ';\n']];
    varList = [varList, ' ', lhs, ' ', rhs];

end


function writeToFile(currentPath, functionContents, varList)

    filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'getFlux.m');

    fid = fopen(filepath, 'w'); 
    if fid == -1
        error(['File not found: ', filepath, ...
            newline, ...
            'Please make sure the folder ', fullfile(currentPath{1}, '..', '..', 'temp'), ...
            ' exists.']) ; 
    end

    fprintf(fid, 'function [v, E1] = getFlux()\n');
    fprintf(fid, ['clearvars', varList, '\n']);
    fprintf(fid, ['syms', varList, '\n']);
    fprintf(fid, functionContents);

    fclose(fid);

end


function [fluxEq, freeEnz]  = runFluxEq(currentPath)
    clear getFlux
    folderPath = fullfile(currentPath{1}, '..', '..', 'temp');
    addpath(folderPath);    
    [fluxEq, freeEnz] = getFlux();
    rmpath(folderPath);

end


function [v, Kpos, Kneg] = parseAllostericFunction(ensemble, rxnI, fluxEq, freeEnz)
    Kpos = {};
    Kneg = {};
    
    syms L_0 n_sub Q v 
    Q = L_0*freeEnz^n_sub;
    
    nPosEffectors = numel(ensemble.populations.models(1).rxnParams(rxnI).KposEff);    
    if nPosEffectors > 0
        effectorList = ensemble.posEffectors{1,1}{rxnI};
        [Kpos, effPosFraction] = parseEffectors(effectorList, nPosEffectors);
        Q = Q * (1 / (1 + str2sym(effPosFraction)))^n_sub;
    end
    
    
    nNegEffetors = numel(ensemble.populations.models(1).rxnParams(rxnI).KnegEff);
    if nNegEffetors > 0
        effectorList = ensemble.negEffectors{1,1}{rxnI};
        [Kneg, effNegFraction] = parseEffectors(effectorList, nNegEffetors);
        Q = Q * (1 + str2sym(effNegFraction))^n_sub;
    end
    
    %syms L eR n Q v 
    %eR = freeEnz;
    %Q = L*freeEnz^n * ((1 + str2sym(effNegFraction))/(1 + str2sym(effPosFraction)))^n;
    v = n_sub*fluxEq / (1 + Q);
    v = simplify(v);
    
end 
    
    
function [Keff, effFraction] = parseEffectors(effectorList, nEffectors)
    
    Keff = {};
    effFraction = '';
    
    for effectorI=1:nEffectors
    
        effector = effectorList{effectorI};
        Keff{effectorI} = ['Keff_', effector];
        effFraction = [effFraction, ' + (', effector, ' / ', Keff{effectorI}, ')'];
    
    end
    
    effFraction = effFraction(2:end);
    
end


function [SBmodel, reactionRate] = addReactionToSBModel(SBmodel, ensemble, fluxEq, modelI, rxnI, Kpos, Kneg, catalyticFxn)
    
    rxnName = ensemble.rxns{rxnI};
    rxnString = generateRxnString(ensemble, rxnI);
    r = addreaction (SBmodel, rxnString);
    rename(r, rxnName);

    SBmodel = addRateConsts(SBmodel, ensemble, modelI, rxnI);
    if isfile(catalyticFxn)
        SBmodel = addAllostericParams(SBmodel, ensemble, modelI, rxnI, Kpos, Kneg);
    end

    reactionRate = regexprep(char(fluxEq), 'k(\d+)', ['k$1_', rxnName]);
    reactionRate = regexprep(reactionRate, 'n_sub', ['n_sub_', rxnName]);
    reactionRate = regexprep(reactionRate, 'L_0', ['L_0_', rxnName]);
    reactionRate = regexprep(reactionRate, 'Keff_(\w+)', ['Keff_$1_', rxnName]);
    SBmodel.Reactions(rxnI).ReactionRate = reactionRate;

end


function rxnString = generateRxnString(ensemble, rxnI)

    substrates = ensemble.mets(find(ensemble.S(:, rxnI)<0));
    subsString = '';
    for metI=1:numel(substrates)
        subsString = [subsString, ' + ',substrates{metI}];
    end
    subsString = subsString(4:end);
    
    products = ensemble.mets(find(ensemble.S(:, rxnI)>0));
    prodString = '';
    for metI=1:numel(products)
        prodString = [prodString, ' + ', products{metI}];
    end
    prodString = prodString(4:end);
    
    rxnString = [subsString, ' <-> ', prodString];

end


function SBmodel = addRefConc(SBmodel, ensemble, modelI)

    for metI=1:numel(ensemble.mets)
        metName = ensemble.mets{metI};
        metConcRef = ensemble.populations.models(modelI).metConcRef(metI);
        param = addparameter(SBmodel, [metName, '_ref'], 'Value', metConcRef * 10^3);
        set(param, 'Units', 'millimol/L');
    end
    
end


function SBmodel = addRateConsts(SBmodel, ensemble, modelI, rxnI)
	
    rateConstNames = getRateConstNames(ensemble, rxnI, '');
    rateConstNames = strsplit(rateConstNames(2:end), ' ');
    
    for rateI=1:numel(rateConstNames)
        rateName = [rateConstNames{rateI}, '_', ensemble.rxns{rxnI}];
        rateVal = ensemble.populations.models(modelI).rxnParams(rxnI).kineticParams(rateI);
        addparameter(SBmodel, rateName, 'Value', rateVal);
    end
    
end


function SBmodel = addAllostericParams(SBmodel, ensemble, modelI, rxnI, Kpos, Kneg)
    
    rxnName = ensemble.rxns{rxnI};
    paramName = ['L_0_', rxnName];
    paramValue = ensemble.populations.models(modelI).rxnParams(rxnI).L;
    addparameter(SBmodel, paramName, 'Value', paramValue);
    
    paramName = ['n_sub_', rxnName];
    paramValue = ensemble.subunits{1}(rxnI);
    addparameter(SBmodel, paramName, 'Value', paramValue);
    
    for KposI=1:numel(Kpos)
        paramName = [Kpos{KposI}, '_', rxnName];
        paramValue = ensemble.populations.models(modelI).rxnParams(rxnI).KposEff(KposI);
        addparameter(SBmodel, paramName, 'Value', paramValue);
    end
    
    for KnegI=1:numel(Kneg)
        paramName = [Kneg{KnegI}, '_', rxnName];
        paramValue = ensemble.populations.models(modelI).rxnParams(rxnI).KnegEff(KnegI);
        addparameter(SBmodel, paramName, 'Value', paramValue);
    end

end


function SBmodel = setSpeciesInitialConcentrations(SBmodel, ensemble)
    
    for metI=1:numel(ensemble.mets)
        met = sbioselect(SBmodel, 'Name', ensemble.mets{metI});
        set(met, 'InitialAmount', 1);
    end
end


function SBmodel = addODEsToSBModel(SBmodel, ensemble,fluxVector)
    
    for metI=1:numel(ensemble.metsBalanced)

        ode = 0;
        metBalI = ensemble.metsBalanced(metI);
        metName = ensemble.mets{metBalI};
        
        consumed = find(ensemble.S(metBalI,:) == -1);
        produced = find(ensemble.S(metBalI,:) == 1);
        
        for consI=1:numel(consumed)
            ode = ode - fluxVector{consumed(consI)};
        end
        
        for prodI=1:numel(produced)
            ode = ode + fluxVector{produced(prodI)};
        end
        
        ode = ode * str2sym(['1 / ', metName, '_ref']);
        %ode = simplify(ode);
        
        rule = addrule(SBmodel, [metName, '= ', char(ode)]);
        %rule = addrule(SBmodel, 'm_m_g6p_c = - 1/m_m_g6p_c * (k05_r_R_FBA1*m_m_g6p_c-k04_r_R_FBA1*m_m_f6p_c)');
        set(rule, 'RuleType', 'rate');
    end
    
    %for metI=1:numel(ensemble.metsBalanced)
    %    
    %    ruleA = addrule(SBModel, 'a = - 1/a_ref * (k1_r1*a-k2_r1*b)');
    %    set(ruleA, 'RuleType', 'rate');
    %end
end


function SBmodel = setConstSpecies(SBmodel, ensemble)
    constMets = setdiff(1:numel(ensemble.mets), ensemble.metsBalanced);
    
    for metI=1:numel(constMets)
        met = sbioselect(SBmodel, 'Name', ensemble.mets{constMets(metI)});
        set(met, 'Constant', true, 'BoundaryCondition', true);
    end
end


function postProcessSBMLFile(SBMLFile, ensemble)
    
    SBMLcontents = '';
    
    if isfile(SBMLFile)
        fileIn = fopen(SBMLFile);
        line = fgets(fileIn);
        SBMLcontents = [SBMLcontents, line];
        
        while ~contains(line, '</sbml>')
            if contains(line, '<reaction')
                
                res = regexp(line, 'name="(\S+)"', 'tokens');
                rxnName = res{1};
                
                rxnInd = find(ismember(ensemble.rxns, rxnName)==1);
                line = fgets(fileIn);
                SBMLcontents = [SBMLcontents, line];
                
                while ~contains(line, '</listOfProducts>')
                    line = fgets(fileIn);
                    SBMLcontents = [SBMLcontents, line];
                end
                
                SBMLcontents = addListOfModifiers(SBMLcontents, ensemble, rxnInd);
                
            end
            line = fgets(fileIn);
            SBMLcontents = [SBMLcontents, line];
        end
        
        fclose(fileIn);
    end
    
    SBMLcontents = substituteIds(SBMLcontents);
    
    fileOut = fopen(SBMLFile, 'w');
    fwrite(fileOut, SBMLcontents);
    fclose(fileOut);
    
end


function SBMLcontents = addListOfModifiers(SBMLcontents, ensemble, rxnInd)
    
    nInhibitors = numel(ensemble.inhibitors{1}{rxnInd});
    nActivators = numel(ensemble.activators{1}{rxnInd});
    nPosEff = numel(ensemble.posEffectors{1}{rxnInd});
    nNegEff = numel(ensemble.negEffectors{1}{rxnInd});
    
    if (nInhibitors > 0) || (nActivators > 0) || (nPosEff > 0) || (nNegEff > 0) 
        SBMLcontents = [SBMLcontents, ['        <listOfModifiers>', newline]];
    else
        return
    end     
        
    if nInhibitors > 0
        for inhibI=1:nInhibitors
            SBMLcontents = [SBMLcontents, ['          <modifierSpeciesReference species="', ensemble.inhibitors{1}{rxnInd}{inhibI}, '"/>', newline]];
        end
    end    
    
    if nActivators > 0
        for activI=1:nActivators
            SBMLcontents = [SBMLcontents, ['          <modifierSpeciesReference species="', ensemble.activators{1}{rxnInd}{activI}, '"/>', newline]];
        end
    end    
    
    if nPosEff > 0
        for posEffI=1:nPosEff
            SBMLcontents = [SBMLcontents, ['          <modifierSpeciesReference species="', ensemble.posEffectors{1}{rxnInd}{posEffI}, '"/>', newline]];
        end
    end
    
    if nNegEff > 0
        for negEffI=1:nNegEff
            SBMLcontents = [SBMLcontents, ['          <modifierSpeciesReference species="', ensemble.negEffectors{1}{rxnInd}{negEffI}, '"/>', newline]];
        end
    end
    
    SBMLcontents = [SBMLcontents, ['        </listOfModifiers>', newline]];
    
end


function SBMLcontents = substituteIds(SBMLcontents)

    idNameMap = regexp(SBMLcontents, 'id="(?<id>\S+)" name="(?<name>\S+)"', 'names');
    
    for idNameI=1:numel(idNameMap)
        SBMLcontents = strrep(SBMLcontents, idNameMap(idNameI).id, idNameMap(idNameI).name);
    end
    
end