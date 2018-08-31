function K = getKernelParam(theta,method,weights)
%--------------------------------------------------------------------------
% Computes global parameter of the dirichlet kernel
%
% Inputs:       theta (Dirichlet variates)
%
% Outputs:      K (kernel param)
%--------------------- Pedro Saa 2016 -------------------------------------
switch method
    case 1                      % Weighted total variance
        D   = size(theta,2);
        Var = zeros(D);
        for ix = 1:D
            for jx = 1:D
                if (ix~=jx)
                    Var(ix,jx) = var(log(theta(:,ix)./theta(:,jx)),weights);
                end
            end
        end
        K = 2*D/sum(Var(:));
    case 2					     % Weighted generalized variance
        D   = size(theta,2);
        Var = zeros(D);
        for ix = 1:D
            for jx = 1:D
                if (ix~=jx)
                    Var(ix,jx) = var(log(theta(:,ix)./theta(:,jx)),weights);
                end
            end
        end
        K = 1/det(Var);
    case 3						 % Pseudo concentration parameter of the MLE estimate
        D        = size(theta,2);
        options  = optimset('TolX',1e-10,'TolFun',1e-8,'MaxIter',1e4,'MaxFunEvals',1e4,'Display','off');
        alpha0   = ones(1,D);
        alphaMLE = exp(fminsearch(@probDirichlet,log(alpha0),options,theta,3,weights'));                        % Compute weighted-MLE for alpha        
        K        = sum(alphaMLE)-D;
end