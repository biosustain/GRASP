function Rev = generalRevSampling(Omega,DGr_RT,numSamples,thinning,burn_in)
% General-use Hit-And-Run for sampling microscopic reversibilities
% 
% ------------------ Pedro Saa 2018 ---------------------------------------
%
% Check input parameters
if (nargin<2)
    disp('Not enough input parameters'); return;
elseif (nargin<3)
    numSamples = 1;
    thinning   = 1e1;
    burn_in    = 1e3;
elseif (nargin<4)
    thinning   = 1e1;
    burn_in    = 1e3;
elseif (nargin<5)
    burn_in    = 1e3;    
end

% General parameters
uTol = 1e-6;

% Define projection parameters
A  = Omega;
b  = DGr_RT;
xcurr = pinv(A)*b;

% General solutions is of the form x = N*alpha + xp. Let us first calculate
% a basis for null(A)
N = null(A);

% Find enclosing box for the extended space
Aeq = A;
beq = b;
lb  = min(DGr_RT)*ones(size(Omega,2),1);
ub  = zeros(size(Omega,2),1);
f   = zeros(size(N,1),1);
x_bnd = zeros(size(N,1),2);
for ix = 1:numel(f)
    
    % Solve min problem
    f(ix) = 1;
    [~,fval] = linprog(f,[],[],Aeq,beq,lb,ub);
    x_bnd(ix,1) = fval;
    
    % Solve max problem
    f(ix) = -1;
    [~,fval] = linprog(f,[],[],Aeq,beq,lb,ub);
    x_bnd(ix,2) = -fval;
    
    % Reset objective
    f(ix) = 0;
end

% Figure out the true max & min step sizes
counter = 0;
xcount  = 0;
Rev     = zeros(size(N,1),numSamples);
while counter < burn_in + numSamples*thinning
    
    % Advance one step
    counter = counter+1;
    
    % Sample direction
    udir = N(:,randi(size(N,2)));
    
    % Figure out where to move
    posDir      = bsxfun(@gt,udir,uTol);
    negDir      = bsxfun(@lt,udir,-uTol);
    distUb      = bsxfun(@plus,-xcurr,x_bnd(:,2));
    distLb      = bsxfun(@minus,xcurr,x_bnd(:,1));
    posStepTemp = bsxfun(@rdivide,distUb,udir);
    negStepTemp = bsxfun(@rdivide,-distLb,udir);
    
    % Figure out step
    posStep = min([posStepTemp(posDir);negStepTemp(negDir)],[],1);
    negStep = max([posStepTemp(negDir);negStepTemp(posDir)],[],1);
    
    % Generate new sample vector
    lambda = negStep + rand(1)*(posStep-negStep);
    xcurr  = xcurr + lambda*udir;
    
    % Save new point
    if (burn_in<=counter) && ~mod(counter-burn_in,thinning)
        xcount = xcount + 1;
        Rev(:,xcount) = xcurr;
    end
end
Rev = exp(Rev);       % Return actual reversibilities (not log)
