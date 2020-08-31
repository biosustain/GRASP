function [vMean,vStd] = computeRobustFluxes(Sflux,xMean,xStd)
% Compute robust fluxes for all reactions by employing consistency
% analysis.
%
% If enough fluxes are known, the remaining are determined in a
% deterministic manner by setting
%
% .. math::
%
%       S_{known} * v_{known} = - S_{unknown} * v_{unknown}
%
% and solving for :math:`v_{unknown}`. 
%
% Otherwise, [TODO : need Pedro's help xD]
%
%
% USAGE:
%
%    [vMean, vStd] = computeRobustFluxes(S, xMean, xStd)
%
% INPUT:
%    Sflux (int matrix):      stoichiometric matrix used for flux calculations
%    xMean (double vector):   mean value for known reactions fluxes
%    xStd (double vector): 	  standard deviations for  known reactions fluxes
%
% OUTPUT:
%    vMean (double matrix):   range of feasible Gibbs energies
%    vStd (double matrix):    range of feasible metabolite concentrations
%    vrng (double matrix):    range of feasible reaction fluxes
%
% .. Authors:
%       - Pedro Saa     2016 original code
%       - Marta Matos   2018 added assertion

% Determine measured fluxes and decompose stoichiometric matrix
idxMeas = find(xMean~=0);
idxUnkn = find(xMean==0);
Sm = Sflux(:,idxMeas);
Sc = Sflux(:,idxUnkn);

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
    
    assert(size(vMean(~vMean), 1) == size(xMean(idxMeas), 1), 'size(vMean(~vMean), 1) ~= size(xMean(idxMeas), 1), most likely some met that should be balanced was set as not balanced or vice versa. Check the met sheet. Also make sure the reactions in measRates are part of the stoichiometric matrix.');    
    
    vMean(~vMean)  = xMean(idxMeas);                                        
    vStd(idxUnkn)  = diag(pinv(Sc)*Sm*Dm*Sm'*pinv(Sc)');
    vStd(~vStd)    = diag(Dm);
    
    % Else, perform gross error analysis
else
    errX = Rred*xMean(idxMeas);    % Compute covariance matrix
    Derr = Rred*Dm*Rred';
    herr = errX'*inv(Derr)*errX;   % Compute statistic and check consistency
    p    = chi2cdf(herr,rank(Rred));
    if (p<5e-2)
        disp('Model is consistent')
    end
    
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
