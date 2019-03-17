% Defining Extreme Pathways
% Input is the elementary matrix
% Output is the Transpose of the flux vector for the extreme pathways

clear cell


S = [-1 1 0 0 0 0 0; -1 0 1 0 0 0 0; 0 -1 0 1 0 0 0; 0 0 -1 1 0 0 0;... 
  0 0 0 -1 1 0 0; 0 0 0 0 -1 1 0; 0 0 0 0 -1 0 1; 1 0 0 0 0 -1 0; 1 0 0 0 0 0 -1];
S = [S; -1.*S];

for xi = 1:size(S, 1)
    s = S(xi,:);
    Q = S == -1.*s;
    q = sum(Q, 2);
    S(q == size(S, 2),:) = [];
    if(xi+1>size(S,1))
        break
    end
end
    
    

I = eye(size(S,1));

A = [I S];

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
B = cell(size(A,1),1);
for o = 1:size(S,2)
    for q = 1:size(A,1)
        if A(q,o) == 0
            B{q} = [B{q}, o];
        end
    end
end

C = [];

for cell = 1:size(B,1)
    for comp = 1:size(B,1)
        if all(ismember(B{cell}, B{comp})) && cell ~= comp
            C = [C, cell];
            B(cell) = {randg};
        end
    end
end

A(C,:) = [];

A = A(:,[1:size(S,1)]);

A'