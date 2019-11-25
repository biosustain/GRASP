function [v,E1,E2] = r_R_GTHPi1(X,K)
% Metabolites definition 
A = X(1);
B = X(2);
C = X(3);
P = X(4);
% Parameters definition K 
k01 = K(1);
k02 = K(2);
k03 = K(3);
k04 = K(4);
k05 = K(5);
k06 = K(6);
k07 = K(7);
k08 = K(8);
k09 = K(9);
k10 = K(10);
%  Numerator terms
E1 = k02*k09*k04*k06+k02*k09*k04*k07+k02*k09*k07*k05.*C+k02*k04*k06*k08+k09*k07*k05.*C*k03.*B;
E2 = k01.*A*k04*k09*k06+k01.*A*k04*k09*k07+k01.*A*k09*k07*k05.*C+k01.*A*k04*k06*k08+k04*k06*k08*k10.*P;
E3 = k03.*B*k06*k01.*A*k09+k03.*B*k01.*A*k09*k07+k06*k08*k10.*P*k02+k03.*B*k06*k01.*A*k08+k03.*B*k06*k08*k10.*P;
E4 = k05.*C*k03.*B*k01.*A*k09+k08*k10.*P*k02*k04+k05.*C*k08*k10.*P*k02+k05.*C*k08*k03.*B*k01.*A+k05.*C*k08*k03.*B*k10.*P;
E5 = k10.*P*k02*k04*k06+k10.*P*k07*k02*k04+k10.*P*k07*k02*k05.*C+k07*k05.*C*k03.*B*k01.*A+k10.*P*k07*k05.*C*k03.*B;
% Denominator terms
D = E1+E2+E3+E4+E5;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
% Reaction rate 
v = +k09.*E5-k10.*P.*E1;
