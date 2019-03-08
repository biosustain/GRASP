function [v,E1,E2] = r_DPGM1(X,K)
% Metabolites definition 
A = X(1,:);
P = X(2,:);
Q = X(3,:);
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
E1 = k02*k07*k04+k02*k07*k05+k02*k04*k06.*P+k07*k05*k03;
E2 = k01.*A*k04*k07+k01.*A*k07*k05+k01.*A*k04*k06.*P+k04*k06.*P*k08.*Q;
E3 = k03*k01.*A*k07+k06.*P*k08.*Q*k02+k03*k06.*P*k01.*A+k03*k06.*P*k08.*Q;
E4 = k08.*Q*k02*k04+k08.*Q*k05*k02+k05*k03*k01.*A+k08.*Q*k05*k03;
% Denominator terms
D = E1+E2+E3+E4;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
% Reaction rate 
v = +k05.*E3-k06.*P.*E4;