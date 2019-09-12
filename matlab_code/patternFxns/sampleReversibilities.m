function randomRev = sampleReversibilities(revMatrix)
%--------------------------------------------------------------------------
% Sample uniformly reversibilities using a Dir(1) distribution
%
% Inputs:     (revMatrix) reversibility matrix with all the cycles in the
%                         reaction pattern
%
% Outputs:    (randomRev) feasible random reversibilities from the pattern
%-----------------------Pedro Saa 2016-------------------------------------
% 1. Solve sequential mechanism
if (size(revMatrix,1)==1)
    tempRev                 = -log(1-rand(sum(revMatrix),1));
    randomRev(revMatrix~=0) = tempRev/sum(tempRev);
    randomRev               = randomRev(:);

% 2. Solve branched mechanism
else

    % We define auxiliar variables to sample the reversibility space
    lbRev     = -log(1-rand(n,1));
    lbRev     = lbRev/sum(lbRev);
    randomRev = computeBrancheReversibilities(revMatrix,lbRev);
end