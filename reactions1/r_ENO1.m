function v = r_ENO1(X,Kr,L,n) 
% Parameters definition 
[vR,eR] = r_ENO1Catalytic(X,Kr); 
Q = L*eR.^n; 
% Reaction rate 
v = n*vR./(1 + Q);