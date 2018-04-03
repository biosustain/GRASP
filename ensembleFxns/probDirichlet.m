function prob = probDirichlet(alpha,x,flag,weights)
%--------------------------------------------------------------------------
% Compute Dir probability or weighted log-likelihood
%
% Inputs:       xdata (in rows), alpha , weights  , flag
%
% Outputs:      prob (f or log-likelihood f)
%--------------------- Pedro Saa 2016 -------------------------------------
% Return probability evaluated at xdata
if (flag==1)
    betaFxn = prod(gamma(alpha))/gamma(sum(alpha));                        % Compute Beta fxn
    prob    = prod(x.^(alpha(ones(size(x,1),1),:)-1),2)/betaFxn;
    
% Return weighted -log:likelihodd
elseif (flag==2)
    betaFxn = prod(gamma(alpha))/gamma(sum(alpha));                        % Compute Beta fxn
    prob    = -weights*log(prod(x.^(alpha(ones(size(x,1),1),:)-1),2)/betaFxn);

% Return weighted -log:likelihodd for unconstrained optimization (alpha' ~ log(alpha)
elseif (flag==3)
    betaFxn = prod(gamma(exp(alpha)))/gamma(sum(exp(alpha)));              % Compute Beta fxn
    prob    = -weights*log(prod(x.^(exp(alpha(ones(size(x,1),1),:))-1),2)/betaFxn);

% Return probability evaluated at (xdata,alpha) (vectorized version)
elseif (flag==4)
	betaFxn = prod(gamma(alpha),2)./gamma(sum(alpha,2));
	prob    = prod(x(ones(size(alpha,1),1),:).^(alpha-1),2)./betaFxn;
end