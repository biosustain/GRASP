function v = r_GAPDH1(X,Kr,L,n) 
% Parameters definition 
[vR,eR] = r_GAPDH1Catalytic(X,Kr); 
Q = L*eR.^n; 
% Reaction rate 
v = n*vR./(1 + Q);