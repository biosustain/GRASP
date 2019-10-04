function [fc,grad] = poolConstraintFxn(x,A,b)
% [TODO: Pedro description]
%
%
% USAGE:
%
%    [fc, grad] = poolConstraintFxn(x, A, b)
%
% INPUT:
%    x (double vector):	  MCA results
%    A (double matrix):   model ensemble
%    b (double vector):   categories to be included in the plot
%
% OUTPUT:
%    fc (double matrix):   MCA results
%    grad (double matrix): model ensemble
%
% .. Authors:
%       - Pedro Saa         2016 original code

fc = (A*x(:) - b);
if (nargout > 1)
    grad = A;	
end