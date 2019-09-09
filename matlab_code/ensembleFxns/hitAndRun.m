function hitAndRun(A, G0Lb,G0Ub,metsLb,metsUb,dGLb,dGUb, S)

RT = 8.314*298.15/1e3; 

numSamples  = 1;
thinning    = 1e1;
burn_in     = 1e3;

nMets = size(metsLb, 1);
nRxns = size(G0Lb, 1);


% Generate a random interior point
posQ = find(metsLb > 0);
negQ = find(metsLb < 0);
sense = zeros(1,nRxns*2);
sense(1:nRxns) = '>';
sense((nRxns + 1):nRxns*2) = '<';

lb = -1000*ones(nMets+nRxns,1);
ub = 100*ones(nMets+nRxns,1);

lb(1:nRxns) = G0Lb;
lb((nRxns+1):(nMets+nRxns)) = metsLb;


ub(1:nRxns) = G0Ub; 
ub((nRxns+1):(nMets+nRxns)) = metsUb;



gurobiModel.A = sparse([eye(nRxns), RT*S';...
                        eye(nRxns), RT*S';]);

gurobiModel.rhs = [dGLb;...
                   dGUb];

               
gurobiModel.sense = char(sense);

gurobiModel.lb = lb;
gurobiModel.ub = ub;
gurobiModel.vtype = 'C';

gurobi_write(gurobiModel, 'try.lp');

sol = gurobi(gurobiModel);



xcurr=sol.x;
lb = [G0Lb; metsLb];
ub = [G0Ub; metsUb];

N = null(A);
uTol = 1e-8;
% Initialize sampling parameters
counter = 0;
xcount  = 0;
Rev     = zeros(size(N,1),numSamples);
while counter < burn_in + numSamples*thinning
    
    % Advance one step
    counter = counter+1;
    
    % Draw direction of movement until a move is possible
    flag = true;
    %while flag
        
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
        %flag = (posStep-negStep < uTol) | (posStep < 0) | (negStep > 0);
    %end
    
    % Generate new sample according to the linear step size fxn
  
    lambda = negStepTemp + rand(1)*(posStepTemp-negStepTemp);
    xcurr  = xcurr + lambda.*udir;
    
    % Ensure the point is within the boundaries
    xcurr(bsxfun(@gt,xcurr,ub)) = ub(bsxfun(@gt,xcurr,ub)); 
    xcurr(bsxfun(@lt,xcurr,lb)) = lb(bsxfun(@lt,xcurr,lb));
    
    bsxfun(@gt, [eye(nRxns), RT*S']*xcurr,dGLb);
    bsxfun(@lt, [eye(nRxns), RT*S']*xcurr,dGUb);
    
    
    % Save new point
    if (burn_in<=counter) && ~mod(counter-burn_in,thinning)
        xcount = xcount + 1;
        Rev(:,xcount) = xcurr;
    end
end

end

