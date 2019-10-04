function vector = colVector(vector)
% Transforms a vector into a column-vector.
%
%
% USAGE:
%
%    vector = colVector(vector)
%
% INPUT:
%    vector (vector):   vector or structure array
%
% OUTPUT:
%    vector (vector):	column-vector or structure array
%
% .. Authors:
%       - Pedro Saa     2014 original code

if size(vector,1)<size(vector,2)
    vector = vector';
end