function [y,U_full,U_ilrt] = map2uncons(x)
x         = (x(:))';
numParams = numel(x);
for jx = 1:numParams-1
    U_ilrt(:,jx) = sqrt(jx/(jx+1))*[(1/jx)*ones(jx,1);-1;zeros(numParams-jx-1,1)];
end
U_full = [U_ilrt,ones(numParams,1)/sqrt(numParams)];
y = log(x)*U_ilrt;