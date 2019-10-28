function branchedRev = computeBranchedReversibilities(revMatrix,lbRev)
% Compute branched reversibilities using our modified approach.
%
% [TODO more details? Pedro?]
%
%
% USAGE:
%
%    branchedRev = computeBranchedReversibilities(revMatrix, lbRev)
%
% INPUT:
%    revMatrix (int matrix):	reversibility matrix with all the cycles in the reaction pattern   
%    lbRev (double vector):	    lower bound for the reversibility values
%
% OUTPUT:
%    branchedRev (double vector):	feasible random reversibilities from the pattern
%
% .. Authors:
%       - Pedro Saa     2016 original code

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