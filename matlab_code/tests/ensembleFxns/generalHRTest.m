classdef generalHRTest < matlab.unittest.TestCase

    properties
        currentPath
    end
    
    methods(TestClassSetup)
        function defineCurrentPath(testCase)
            testCase.currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
        end
    end
    
 
    methods (Test)
         function testGeneralHRFluxesNormal(testCase)
            
            seed = 1;
            rng(seed);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs.mat'));
            ensemble = ensemble.ensemble;
            
            nSamples= 100;
            nSteps = 1e2;
            nDiscard= 1000;
            priorType = 'normal';
            
            fluxLB     = ensemble.fluxRanges(:,1);
            fluxUB     = ensemble.fluxRanges(:,2);
            fluxX0     = ensemble.initialTMFAPoint(1:size(fluxLB,1));
            fluxAeq    = ensemble.Sflux;
            fluxPoints = generalHR(fluxAeq,fluxLB,fluxUB,fluxX0,nSamples,nSteps,nDiscard,priorType);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_generalHRTestFluxesNormal.mat'));
            trueRes = trueRes.fluxPoints;
            
            testCase.verifyThat(fluxPoints, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)))
         end
        
         function testGeneralHRFluxesUniform(testCase)
            
            seed = 1;
            rng(seed);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs.mat'));
            ensemble = ensemble.ensemble;
            
            nSamples= 100;
            nSteps = 1e2;
            nDiscard= 1000;
            priorType = 'uniform';
            
            fluxLB     = ensemble.fluxRanges(:,1);
            fluxUB     = ensemble.fluxRanges(:,2);
            fluxX0     = ensemble.initialTMFAPoint(1:size(fluxLB,1));
            fluxAeq    = ensemble.Sflux;
            fluxPoints = generalHR(fluxAeq,fluxLB,fluxUB,fluxX0,nSamples,nSteps,nDiscard,priorType);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_generalHRTestFluxesUniform.mat'));
            trueRes = trueRes.fluxPoints;
            
            testCase.verifyThat(fluxPoints, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)))
        end
        
        function testGeneralHRFGNormal(testCase)
            
            seed = 1;
            rng(seed);
           
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs.mat'));
            ensemble = ensemble.ensemble;
            
            nSamples= 100;
            nSteps = 1e2;
            nDiscard= 1000;
            priorType = 'normal';
            
            RT           = 8.314*298.15/1e3;  % [kJ/mol]
            [m, n]       = size(ensemble.Sthermo);
            thermoLB     = [ensemble.gibbsRanges(ensemble.idxNotExch,1);ensemble.DGfStdRange(:,1);ensemble.lnMetRanges(:,1)];
            thermoUB     = [ensemble.gibbsRanges(ensemble.idxNotExch,2);ensemble.DGfStdRange(:,2);ensemble.lnMetRanges(:,2)];
            thermoAeq    = [eye(size(ensemble.Sthermo,2)),-ensemble.Sthermo',-RT*ensemble.Sthermo'];
            thermoX0     = ensemble.initialTMFAPoint(n+1:end);
            thermoPoints = generalHR(thermoAeq,thermoLB,thermoUB,thermoX0,nSamples,nSteps,nDiscard,priorType);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_generalHRTestDGNormal.mat'));
            trueRes = trueRes.thermoPoints;
            
            testCase.verifyThat(thermoPoints, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)))
        end
        
        function testGeneralHRFGUniform(testCase)
            
            seed = 1;
            rng(seed);
            
            ensemble = load(fullfile(testCase.currentPath{1}, 'testFiles', 'initializedEnsemble_toy_model1_sampleFluxGibbs.mat'));
            ensemble = ensemble.ensemble;
           
            nSamples= 100;
            nSteps = 1e2;
            nDiscard= 1000;
            priorType = 'uniform';            
            
            RT           = 8.314*298.15/1e3;  % [kJ/mol]
            [m, n]       = size(ensemble.Sthermo);
            thermoLB     = [ensemble.gibbsRanges(ensemble.idxNotExch,1);ensemble.DGfStdRange(:,1);ensemble.lnMetRanges(:,1)];
            thermoUB     = [ensemble.gibbsRanges(ensemble.idxNotExch,2);ensemble.DGfStdRange(:,2);ensemble.lnMetRanges(:,2)];
            thermoAeq    = [eye(size(ensemble.Sthermo,2)),-ensemble.Sthermo',-RT*ensemble.Sthermo'];
            thermoX0     = ensemble.initialTMFAPoint(n+1:end);
            thermoPoints = generalHR(thermoAeq,thermoLB,thermoUB,thermoX0,nSamples,nSteps,nDiscard,priorType);

            trueRes = load(fullfile(testCase.currentPath{1}, 'testFiles', 'trueRes_generalHRTestDGUniform.mat'));
            trueRes = trueRes.thermoPoints;
            
            testCase.verifyThat(thermoPoints, matlab.unittest.constraints.IsEqualTo(trueRes, ...
                'Within', matlab.unittest.constraints.RelativeTolerance(1e-12)))
        end
    end
end