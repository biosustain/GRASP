function vector = colVector(vector)
%--------------------------------------------------------------------------
% Transforms a vector into a column-vector
%
% Inputs:       (vector)  vector or structure array
%
% Outputs:      (vector)  column-vector or structure array
%-----------------------Pedro Saa 2014-------------------------------------
if size(vector,1)<size(vector,2)
    vector = vector';
end