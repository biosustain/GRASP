function log_p = logp_ilrtMKG(z,alpha,U_full)

% Get number of parameters
numParams = numel(z)+1;
        
% Formulate determinant of the transformation's Jacobian for the new density (we map z onto the full space)
detJ = det(U_full.*exp(U_full*repmat((log( map2simplex(z) )*U_full)',1,numParams)));
        
% Formulate transformed log-prior density (we map z onto the full space)
%               log( exp(             sum( )       )         )                      +  log( det_term ) -     log ( betaFxn )
log_p = (alpha-1)*(U_full*(log( map2simplex(z) )*U_full)') + log(detJ) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));