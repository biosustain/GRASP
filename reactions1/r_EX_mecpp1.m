function [v,E1,E2] = r_EX_mecpp1(X,K)
% Metabolites definition 
A = X(1,:);
P = X(2,:);
% Parameters definition K 
k01 = K(1);
k02 = K(2);
k03 = K(3);
k04 = K(4);
k05 = K(5);
k06 = K(6);
%  Numerator terms
E1 = k02*k05+k02*k04+k05*k03;
E2 = k01.*A*k05+k01.*A*k04+k04*k06.*P;
E3 = k06.*P*k02+k03*k01.*A+k06.*P*k03;
% Denominator terms
D = E1+E2+E3;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
% Reaction rate 
v = +k05.*E3-k06.*P.*E1;