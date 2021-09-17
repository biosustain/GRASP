classdef simulateEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
        relTol = 1e-8;
        absTol = 1e-8;
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
 
    methods(TestMethodTeardown)
        function removeReactionsFolder(testCase)           

            reactionsFolder = fullfile(testCase.currentPath{1}, '..', '..', '..', 'reactions', 'toy_model1_random2_1');
            
            if exist(reactionsFolder, 'dir')
                rmdir(reactionsFolder, 's');
            end
        end
    end
       
    
    methods (Test)
       function testSimulateEnsemble1(testCase)
                       
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            ensemble.freeVars{end+1} = 'r_r13';
            interruptTime = 40;

            metsAbsOrRel = 'rel';
            
            % Change initial conditions here if you want
            enzymesIC = {{'r2', 1.5}};       % Always relative concentrations for enzymes
            metsIC = {{'m9', 2}};        % Absolute concentrations must be given in mmol/L

            % Specifiy the time of simulation (probably in hours)
            finalTime = 1;
            numModels = 5;
            numCores = 2;
            
            simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime,numModels, numCores);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSimulateEnsemble1'));
            trueRes = trueRes.simulationRes;
                   
            testCase.verifyThat(simulationRes, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
       end
        
        function testSimulateEnsembleAbs(testCase)
                        
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            
            ensemble.freeVars{end+1} = 'r_r13';
            
            interruptTime = 40;

            metsAbsOrRel = 'abs';
            
            % Change initial conditions here if you want
            enzymesIC = {{'r2', 10}};                           % Always relative concentrations for enzymes
            metsIC = {{'m9', 2* 10^-3 }, {'m11', 2 * 10^-1}};        % Absolute concentrations must be given in mmol/L

            % Specifiy the time of simulation (probably in hours)
            finalTime = 1;
            numModels = 5;
            numCores = 2;
            
            simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime, numModels, numCores);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSimulateEnsembleAbs'));
            trueRes = trueRes.simulationRes;
                   

            testCase.verifyThat(simulationRes, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(testCase.relTol) | matlab.unittest.constraints.AbsoluteTolerance(testCase.absTol)));  
        end
    end
end
