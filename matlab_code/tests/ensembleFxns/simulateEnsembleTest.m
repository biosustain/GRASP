classdef simulateEnsembleTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
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
            
            interruptTime = 40;

            % Get default initial conditions, all ones
            freeVars = numel(ensemble.freeVars);
            xopt = ones(freeVars,1);
            ix_mets = 1:numel(ensemble.metsActive);
            ix_enz = ix_mets(end)+1:freeVars;
            metsIC = xopt(ix_mets);
            enzymesIC = xopt(ix_enz);

            metsAbsOrRel = 'rel';
            
            % Change initial conditions here if you want
            enzymesIC(2) = 1.5;
            metsIC(5) = 2;

            % Specifiy the time of simulation (probably in hours)
            finalTime = 1;
            
            simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSimulateEnsemble1'));
            trueRes = trueRes.simulationRes;
                   
            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(simulationRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
       end
        
        function testSimulateEnsembleAbs(testCase)
                        
            % To make sure reaction files are created properly
            xlsxFile = fullfile(testCase.currentPath{1}, 'testFiles', 'toy_model1_random2');
            loadEnsembleStructure(xlsxFile);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'final_ensemble_toy_model1_random2.mat'));
            ensemble = ensemble.ensemble;
            
            interruptTime = 40;

            % Get default initial conditions, all ones
            freeVars = numel(ensemble.freeVars);
            xopt = ones(freeVars,1);
            ix_mets = 1:numel(ensemble.metsActive);
            ix_enz = ix_mets(end)+1:freeVars;
            metsIC = xopt(ix_mets);
            enzymesIC = xopt(ix_enz);

            metsAbsOrRel = 'abs';
            
            % Change initial conditions here if you want
            enzymesIC(2) = 1.5;
            metsIC(5) = 10^-3;
            metsIC(7) = 2 * 10^-1;

            % Specifiy the time of simulation (probably in hours)
            finalTime = 1;
            
            simulationRes = simulateEnsemble(ensemble, finalTime, enzymesIC, metsIC, metsAbsOrRel, interruptTime);
            
            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueResSimulateEnsembleAbs'));
            trueRes = trueRes.simulationRes;
                   

            testCase.verifyThat(trueRes, matlab.unittest.constraints.IsEqualTo(simulationRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-4)));
        end
    end
end
