function Rev = generalRevSampling_new(Omega,DGr_RT,numSamples,thinning,burn_in,xcurr)
% General-purpose reversility sampling
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
    xcurr       = [];
elseif (nargin<4)
    thinning    = 1e1;
    burn_in     = 1e3;
    xcurr       = [];
elseif (nargin<5)
    burn_in     = 1e3;
    xcurr       = [];
elseif (nargin<6)
    xcurr       = [];
end

% Reversibilit tolerance
uTol = 1e-8;

% Find enclosing box for the sampling space. Given the features of the
% space, we can calcultae the hyperbox without running optimizations
LB = Omega.*(repmat(DGr_RT,1,size(Omega,2)));
LB(~LB) = -inf;
lb      = max(LB,[],1)';
lb(isinf(lb)) = 0;
ub = zeros(size(lb));

% General solutions is of the form x = N*alpha + xp. Calculate an orthogonal
% basis for null(A)
N = null(Omega);

% Generate a random interior point
if isempty(xcurr)
    Pn = N*N';
    xp = pinv(Omega)*DGr_RT;
    xcurr = lb/2;
    while true
        xcurr = xp + Pn*xcurr;
        if any(xcurr>uTol) || any(xcurr-lb<-uTol)
            xcurr(xcurr>uTol)     = lb(xcurr>uTol)/2;
            xcurr(xcurr-lb<-uTol) = lb(xcurr-lb<-uTol)/2;
        else
            xcurr(abs(xcurr)<uTol) = 0;
            xcurr(abs(xcurr-lb)<uTol) = lb(abs(xcurr-lb)<uTol);
            break;
        end
    end
end

% Definition of linear step size fxn
stepSizeFxn = @(t0,L,R) L*( 1 + (2*t0-1) - sqrt( (1 + (2*t0-1))^2 - 4*(2*t0-1)*R ) )/(2*(2*t0-1)) - t0;

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
        
        % Sample direction
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
    
    % Generate new sample according to the linear step size fxn
    lambda = stepSizeFxn(-negStep,(posStep-negStep),rand(1));
    xcurr  = xcurr + lambda*udir;
    
    % Ensure the point is within the boundaries
    xcurr(bsxfun(@gt,xcurr,ub)) = ub(bsxfun(@gt,xcurr,ub)); 
    xcurr(bsxfun(@lt,xcurr,lb)) = lb(bsxfun(@lt,xcurr,lb));
    
    % Save new point
    if (burn_in<=counter) && ~mod(counter-burn_in,thinning)
        xcount = xcount + 1;
        Rev(:,xcount) = xcurr;
    end
end
Rev = exp(Rev);       % Return actual reversibilities (not log)