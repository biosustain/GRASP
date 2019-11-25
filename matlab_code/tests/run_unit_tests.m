clearvars
import matlab.unittest.TestSuite

suiteAnalysisFxns = TestSuite.fromFolder('analysisFxns');
suiteEnsembleFxns = TestSuite.fromFolder('ensembleFxns');
suitePatternFxns = TestSuite.fromFolder('patternFxns');

addpath(fullfile('..',  'analysisFxns'), ...
        fullfile('..',  'ensembleFxns'), ...
        fullfile('..',  'patternFxns'));


resultAnalysisFxns = run(suiteAnalysisFxns);
resultEnsembleFxns = run(suiteEnsembleFxns);
resultPatternFxns = run(suitePatternFxns);


disp([newline, newline, 'Result summary:', newline, newline])
disp('Analysis functions')
table(resultAnalysisFxns)
disp('Ensemble functions')
table(resultEnsembleFxns)
disp('Pattern functions')
table(resultPatternFxns)