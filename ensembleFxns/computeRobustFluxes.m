function [vMean,vStd] = computeRobustFluxes(S,xMean,xStd)
%--------------------------------------------------------------------------
% Compute robust fluxes employing consistency analysis
% Inputs:       Stoichiometric matrix, xMeas, xStd
%
% Outputs:      ensemble basic data structure
%
%------------------------Pedro Saa 2016------------------------------------
% Determine measured fluxes and decompose stoichiometric matrix
idxMeas = find(xMean~=0);
idxUnkn = find(xMean==0);
Sm = S(:,idxMeas);
Sc = S(:,idxUnkn);

% Initialize final fluxes
vMean = zeros(size(xMean));
vStd  = zeros(size(xMean));

% Compute estimate Rred
Dm       = diag(xStd(idxMeas).^2);
Rred     = Sm-Sc*pinv(Sc)*Sm;
singVals = svd(Rred);
Rred     = Rred(abs(singVals)>1e-12,:);

% If the system is fully determined, compute as follows
if isempty(Rred)
    vMean(idxUnkn) = -pinv(Sc)*Sm*xMean(idxMeas);
    assert(size(vMean(~vMean), 1) == size(xMean(idxMeas), 1), 'size(vMean(~vMean), 1) ~= size(xMean(idxMeas), 1), most likely some met that should be balanced was set as not balanced or vice versa. Check the met sheet.');    
    vMean(~vMean)  = xMean(idxMeas);                                        
    vStd(idxUnkn)  = diag(pinv(Sc)*Sm*Dm*Sm'*pinv(Sc)');
    vStd(~vStd)    = diag(Dm);
    
    % Else, perform gross error analysis
else
    errX = Rred*xMean(idxMeas);    % Compute covariance matrix
    Derr = Rred*Dm*Rred';
    herr = errX'*inv(Derr)*errX;   % Compute statistic and check consistency
    p    = chi2cdf(herr,rank(Rred));
    if (p<5e-2); disp('Model is consistent'); end;
    
    % Compute balanced measured flux and covariance matrix
    Rp        = Rred'*inv(Derr)*Rred;
    xMean_adj = (eye(size(Dm))-Dm*Rp)*xMean(idxMeas);
    xStd_adj  = Dm-Dm*Rp*Dm;   
    vMean(idxUnkn) = -pinv(Sc)*Sm*xMean_adj;    % Adjust final fluxes
    vMean(~vMean)  = xMean_adj;
    vStd(idxUnkn)  = diag(pinv(Sc)*Sm*xStd_adj*Sm'*pinv(Sc)');
    vStd(~vStd)    = diag(xStd_adj);
end
vStd = sqrt(vStd);                              % Compute std
