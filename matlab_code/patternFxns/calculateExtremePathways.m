function [extremePath] = calculateExtremePathways(Selem)
% This function constructs the extreme pathways of the elementary flux
% matrix.
%
%  [TODO add ref for base paper Nick]
%
%
% USAGE:
%
%    extremePath = calculateExtremePathways(Selem)
%
% INPUT:
%    Selem (int matrix):	Elementary Flux Matrix  
%
% OUTPUT:
%    extremePath (int matrix):	transpose of the flux vector for the extreme pathways
%
% .. Authors:
%       - Nicholas Cowie	2019 original code



S = Selem';

S = [S; -1*S];

% Removes the reverse of each flux in order to ensure that the flux is only
% going in the forward direction
for xi = 1:size(S, 1)
    s = S(xi,:);
    Q = S == -1.*s;
    q = sum(Q, 2);
    S(q == size(S, 2),:) = [];
    if(xi+1>size(S,1))
        break
    end
end
    
    
% Appending the identity matrix with as many rows as metabolites
I = eye(size(S,1));

A = [I S];

% For each metabolite, ensures that the combination of reactions create a
% steady state solution.
% The row operations required to do this are stored in the identiy matrix
% which will become the basis vectors for the elementary rays of the flux
% cone.
for j = size(S,1)+1:size(S,1)+size(S,2)
    A_temp = A;
    A = A_temp(A_temp(:,j) == 0,:);
    Plist = A_temp(:,j) == 1;
    Plist = double(Plist);
    Nlist = A_temp(:,j) == -1;
    Nlist = double(Nlist);
    for x = 1:size(A_temp,1)
        if Plist(x) == 1
            Plist(x) = x;
        end
    end
    
    for x = 1:size(A_temp,1)
        if Nlist(x) == 1
            Nlist(x) = x;
        end
    end
    
    % Creating a vector of indicies for positive and negative values
    Plist = Plist(Plist~=0);
    Nlist = Nlist(Nlist~=0);
    
    for k = Plist'
        for l = Nlist'
            
            A = [A; A_temp(k,:)+A_temp(l,:)];
            
        end
    end
end


% Removing all rows that have a set of zeros as a subset of another row
cell_dim = size(A,1);

B = cell(cell_dim,1);
for o = 1:size(S,2)
    for q = 1:size(A,1)
        if A(q,o) == 0
            B{q} = [B{q}, o];
        end
    end
end

C = [];

for cells = 1:size(B,1)
    for comp = 1:size(B,1)
        if all(ismember(B{cells}, B{comp})) && cells ~= comp
            C = [C, cells];
            B(cells) = {randg};
        end
    end
end

A(C,:) = [];

A = A(:,[1:size(S,1)]);

extremePath = A';
end
