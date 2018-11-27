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
[m,n]       = size(revMatrix);
ub          = ones(n,1);
beq         = ones(m,1);
try
    options.Display = 'off';
    branchedRev = linprog(-ones(n,1),[],[],revMatrix,beq,lbRev,ub,options);
catch
    options     = optimset('Display','off');
    branchedRev = linprog(-ones(n,1),[],[],revMatrix,beq,lbRev,ub,[],options);
end