clearvars
import matlab.unittest.TestSuite

suitePatternFxns = TestSuite.fromFolder('patternFxns');
suiteEnsembleFxns = TestSuite.fromFolder('ensembleFxns');

addpath(fullfile('..',  'patternFxns'), ...
        fullfile('..',  'ensembleFxns'));


resultPatternFxns = run(suitePatternFxns);
resultEnsembleFxns = run(suiteEnsembleFxns);


disp([newline, newline, 'Result summary:', newline, newline])
disp('Pattern functions')
table(resultPatternFxns)
disp('Ensemble functions')
table(resultEnsembleFxns)