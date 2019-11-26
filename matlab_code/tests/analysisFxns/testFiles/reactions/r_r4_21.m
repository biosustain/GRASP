function [v,E1,E2] = r_r4_21(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
P = X(3,:);
Q = X(4,:);
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
E1 = k02*k09*k04*k06+k02*k09*k04*k07+k02*k09*k07*k05+k02*k04*k06*k08.*P+k09*k07*k05*k03.*B;
E2 = k01.*A*k04*k09*k06+k01.*A*k04*k09*k07+k01.*A*k09*k07*k05+k01.*A*k04*k06*k08.*P+k04*k06*k08.*P*k10.*Q;
E3 = k03.*B*k06*k01.*A*k09+k03.*B*k01.*A*k09*k07+k06*k08.*P*k10.*Q*k02+k03.*B*k06*k01.*A*k08.*P+k03.*B*k06*k08.*P*k10.*Q;
E4 = k05*k03.*B*k01.*A*k09+k08.*P*k10.*Q*k02*k04+k05*k08.*P*k10.*Q*k02+k05*k08.*P*k03.*B*k01.*A+k05*k08.*P*k03.*B*k10.*Q;
E5 = k10.*Q*k02*k04*k06+k10.*Q*k07*k02*k04+k10.*Q*k07*k02*k05+k07*k05*k03.*B*k01.*A+k10.*Q*k07*k05*k03.*B;
% Denominator terms
D = E1+E2+E3+E4+E5;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
% Reaction rate 
v = +k07.*E4-k08.*P.*E5;