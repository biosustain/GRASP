function v = r_FBA1(X,Kr,L,n) 
% Parameters definition 
[vR,eR] = r_FBA1Catalytic(X,Kr); 
Q = L*eR.^n; 
% Reaction rate 
v = n*vR./(1 + Q);