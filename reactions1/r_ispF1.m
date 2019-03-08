function v = r_ispF1(X,posEff,Kr,KposEff,L,n) 
% Parameters definition 
[vR,eR] = r_ispF1Catalytic(X,Kr); 
Q = L*eR.^n; 
KposEff = KposEff(ones(size(posEff,2),1),:); 
Q = Q.*((1 + sum(posEff'./KposEff,2)).^-n)'; 
% Reaction rate 
v = n*vR./(1 + Q);