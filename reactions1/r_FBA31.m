function [v,E1,E2] = r_FBA31(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
P = X(3,:);
% Parameters definition K 
k01 = K(1);
k02 = K(2);
k03 = K(3);
k04 = K(4);
k05 = K(5);
k06 = K(6);
k07 = K(7);
k08 = K(8);
%  Numerator terms
E1 = k02*k07*k04+k02*k07*k05+k02*k04*k06+k07*k05*k03.*B;
E2 = k01.*A*k04*k07+k01.*A*k07*k05+k01.*A*k04*k06+k04*k06*k08.*P;
E3 = k03.*B*k01.*A*k07+k06*k08.*P*k02+k03.*B*k06*k01.*A+k03.*B*k06*k08.*P;
E4 = k08.*P*k02*k04+k08.*P*k05*k02+k05*k03.*B*k01.*A+k08.*P*k05*k03.*B;
% Denominator terms
D = E1+E2+E3+E4;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
% Reaction rate 
v = +k07.*E4-k08.*P.*E1;