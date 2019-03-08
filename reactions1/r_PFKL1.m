function v = r_PFKL1(X,negEff,posEff,Kr,KnegEff,KposEff,L,n) 
% Parameters definition 
[vR,eR] = r_PFKL1Catalytic(X,Kr); 
Q = L*eR.^n; 
KnegEff = KnegEff(ones(size(negEff,2),1),:); 
Q = Q.*((1 + sum(negEff'./KnegEff,2)).^n)'; 
KposEff = KposEff(ones(size(posEff,2),1),:); 
Q = Q.*((1 + sum(posEff'./KposEff,2)).^-n)'; 
% Reaction rate 
v = n*vR./(1 + Q);