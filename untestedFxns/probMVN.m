function prob = probMVN(alpha,x)
%--------------------------------------------------------------------------
% Compute MVN probability or weighted log-likelihood
%
% Inputs:       alpha (in cols) , x (query point)
%
% Outputs:      prob (f or log-likelihood f )
%--------------------- Pedro Saa 2016 -------------------------------------
% Get hyperparameters, mu and sigma
mu    = alpha(:,1)';
sigma = alpha(:,2:end);

% Get scalar mean, and use it to center data
d  = size(x,2);
X0 = x - mu(ones(size(x,1),1),:);

% Compute cholesky decomposition of the cov matrix
[R,err] = cholcov(sigma,0);

% Create array of standardized data, and compute log(sqrt(det(sigma)))
xRinv = X0 / R;
logSqrtDetsigma = sum(log(diag(R)));

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2,2);

% Return probability density
prob = exp(-0.5*quadform - logSqrtDetsigma - d*log(2*pi)/2);