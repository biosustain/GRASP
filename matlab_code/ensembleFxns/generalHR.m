function points = generalHR(Aeq,LB,UB,x0,nSamples,nSteps,nDiscard,priorType)
% General Hit and Run algorithm v2.0 (homogeneous version)
%
% Samples from A*x = 0, with lb <= x <= ub 
% This is a different flavour of HR, which generates though the same
% results.
%
%
% USAGE:
%
%    points = generalHR(Aeq,LB,UB,x0,nSamples,nSteps,nDiscard,priorType)
%
% INPUT:
%    Aeq (double vector):        matrix with constraints
%    LB (double vector):         lower bounds for the variables
%    UB (double vector):         upper bounds for the variables
%    x0 (double vector):         initial point to start the sampling (comes from TMFA)
%    nSamples (double vector):   number of samples to get
%    nSteps (double vector):     number of steps
%    nDiscard (double vector):	 number of samples to discard
%    priorType (double vector):  prior type to use in the sampling, either 'normal' or 'uniform'
%
%
% .. Authors:
%       - Pedro A. Saa   	2020 original code

if (nargin<8)
    priorType  = 'uniform';
end

% Determine the dimension of the feasible space
A     = Aeq;
ix_fixed = (abs(UB-LB)==0);
A(:,ix_fixed) = [];
N     = null(A,'r');                                        % get a near sparse null space basis
N     = bsxfun(@rdivide,N,sqrt(sum(N.^2,1)));               % make unitary
fVars = size(N,2);
if (fVars==0)
    points = repmat(x0,1,nSamples);
    return
end

% Define appropriate prior
if strcmp(priorType,'normal')
    range = UB(~ix_fixed)-LB(~ix_fixed);
    Mu    = LB(~ix_fixed) + .5*range;
    R     = range/4;                                        % 99.9% of the support is within +/- 3 std's
%     logSqrtDetSigma = sum(log(R));
    Ln_priorFxn = @(point) -0.5*sum(((point(~ix_fixed) - Mu)./R).^2); % The rest of term is unnecessary as it is constant; -(logSqrtDetSigma+fVars*log(2*pi)/2);
else
    Ln_priorFxn = @(point) 0;      % Uniform prior    
end

% General Hit And Run Sampler
colA        = size(Aeq,2);
uTol        = 1e-9;
bTol        = 1e-9;
currPoint   = x0;                      % Initial point
Ln_currProb = Ln_priorFxn(currPoint);
points      = zeros(colA,nSamples);
count       = 0;
iter        = 0;

% Initiate sampler
while (count<nSamples)
    iter = iter+1;
    
    while true
        
        % Draw sparse random direction. Only draw random directions in dim(Omega)
        udir = zeros(colA,1);
        udir(~ix_fixed) = N(:,randi(fVars));        % draw random unitary direction from N
        
        % Figure out max distance
        posDir      = (udir>uTol);
        negDir      = (udir<-uTol);
        posStepTemp = (UB-currPoint)./udir;
        negStepTemp = (LB-currPoint)./udir;
        posStep     = min([posStepTemp(posDir);negStepTemp(negDir)]);
        negStep     = max([negStepTemp(posDir);posStepTemp(negDir)]);
        
        % Draw new direction if too close to the boundaries
        Lcord = posStep-negStep;
        flag  = (Lcord < bTol) | (posStep < 0) | (negStep > 0);
        if ~flag; break; end
    end
    
    % Sample step size and perform random jump
    lambda    = negStep + rand(1)*Lcord;
    nextPoint = currPoint + lambda*udir;
    
    % Accept/Reject proposal using MCMC kernel
    Ln_nextProb   = Ln_priorFxn(nextPoint);
    Ln_acceptProb = Ln_nextProb - Ln_currProb;
    if (Ln_acceptProb < log(rand(1)))
        nextPoint = currPoint;
    else
        Ln_currProb = Ln_nextProb;
    end
    
    % Update current point
    currPoint = nextPoint;
    
    % Save next point
    if (nDiscard <= iter) && ~mod(iter-nDiscard,nSteps)
        count = count + 1;
        points(:,count) = nextPoint;
    end
end
