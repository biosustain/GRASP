function branchedRev = computeBranchedReversibilities(revMatrix,lbRev)
%--------------------------------------------------------------------------
% Compute branched reversibilities using our modified approach
%
% Inputs:     (revMatrix) reversibility matrix with all the cycles in the
%                         reaction pattern
%             (lbRev)     branched
%
% Outputs:    (randomRev) feasible random reversibilities from the pattern
%-----------------------Pedro Saa 2016-------------------------------------
% Cast optimization problem and solve branched Rev
options     = optimset('Display','off');
[m,n]       = size(revMatrix);
ub          = ones(n,1);
beq         = ones(m,1);
branchedRev = linprog(-ones(n,1),[],[],revMatrix,beq,lbRev,ub,[],options);
% branchedRev = linprog(-ones(n,1),[],[],revMatrix,beq,lbRev,ub);