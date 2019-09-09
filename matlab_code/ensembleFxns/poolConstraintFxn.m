function [fc,grad] = poolConstraintFxn(x,A,b)
fc = (A*x(:) - b);
if (nargout > 1)
    grad = A;	
end