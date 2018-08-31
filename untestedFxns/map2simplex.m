function [x,U_full,U_ilrt] = map2simplex(y)
y         = (y(:))';
numParams = numel(y)+1;
for jx = 1:numParams-1
    U_ilrt(:,jx) = sqrt(jx/(jx+1))*[(1/jx)*ones(jx,1);-1;zeros(numParams-jx-1,1)];
end
U_full = [U_ilrt,ones(numParams,1)/sqrt(numParams)];
x      = (exp(U_ilrt*y'))'./repmat(sum((exp(U_ilrt*y'))',2),1,numParams);