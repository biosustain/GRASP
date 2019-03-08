function v = r_PGI1(X,Kr,L,n) 
% Parameters definition 
[vR,eR] = r_PGI1Catalytic(X,Kr); 
Q = L*eR.^n; 
% Reaction rate 
v = n*vR./(1 + Q);