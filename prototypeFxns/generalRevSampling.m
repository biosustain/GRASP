function Rev = generalRevSampling(Omega,DGr_RT,numSamples,thinning,burn_in)
% Thermodynamic-based Flux Balance Analysis
% 
% ------------------ Pedro Saa 2018 ---------------------------------------
%
% Check input parameters
if (nargin<2)
    disp('Not enough input parameters'); return;
elseif (nargin<3)
    numSamples  = 1;
    thinning    = 1e1;
    burn_in     = 1e3;
elseif (nargin<4)
    thinning    = 1e1;
    burn_in     = 1e3;
elseif (nargin<5)
    burn_in     = 1e3;
end

% General parameters
uTol = 1e-6;

% Define projection parameters
A  = Omega;
b  = DGr_RT;

% General solutions is of the form x = N*alpha + xp. Calculate a parse
% basis for null(A) that favours movement
N = null(A,'r');
N = bsxfun(@rdivide,N,sqrt(sum(N.^2)));

% Find enclosing box for the sampling space. Given the features of the
% space, we can calcultae the hyperbox without running optimizations
LB = A.*(repmat(b,1,size(A,2)));
LB(~LB) = -inf;
lb      = max(LB,[],1)';
lb(isinf(lb)) = 0;
Aeq = A;
beq = b;
ub  = zeros(size(Omega,2),1);

% Solve 1 optimization with random objective and force it to the interior.
% Use this solution as initial point
alpha = .1;
xcurr = linprog(-rand(size(lb))./abs(lb),[],[],Aeq,beq,(1-alpha)*lb,alpha*lb);

% Initialize sampling parameters
counter = 0;
xcount  = 0;
Rev     = zeros(size(N,1),numSamples);
while counter < burn_in + numSamples*thinning
    
    % Advance one step
    counter = counter+1;
    
    % Draw direction of movement until a move is possible
    flag = true;
    while flag
        
        % Sample direction N*rand(size(N,2),1);%
        udir = N(:,randi(size(N,2)));
        
        % Figure out where to move
        posDir      = bsxfun(@gt,udir,uTol);
        negDir      = bsxfun(@lt,udir,-uTol);
        distUb      = bsxfun(@plus,-xcurr,ub);
        distLb      = bsxfun(@minus,xcurr,lb);
        posStepTemp = bsxfun(@rdivide,distUb,udir);
        negStepTemp = bsxfun(@rdivide,-distLb,udir);
        
        % Figure out step
        posStep = min([posStepTemp(posDir);negStepTemp(negDir)],[],1);
        negStep = max([posStepTemp(negDir);negStepTemp(posDir)],[],1);
        
        % Draw new direction if too close to the boundaries
        flag = (posStep-negStep < uTol) | (posStep < 0) | (negStep > 0);
    end
    
    % Generate new sample according to the appropriate step size fxn
    lambda = negStep + rand(1)*(posStep-negStep);
    xcurr  = xcurr + lambda*udir;
    
    % Save new point
    if (burn_in<=counter) && ~mod(counter-burn_in,thinning)
        xcount = xcount + 1;
        Rev(:,xcount) = xcurr;
    end
end
Rev = exp(Rev);       % Return actual reversibilities (not log)
