function [v,E1,E2] = r_dxr1(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
I = X(3,:);
P = X(4,:);
Q = X(5,:);
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
k11 = K(11);
k12 = K(12);
%  Numerator terms
E1 = k02*k11*k04*k10*k06+k02*k11*k04*k07*k10+k02*k11*k07*k10*k05+k02*k04*k06*k08.*P*k10+k11*k07*k10*k05*k03.*B;
E2 = k01.*A*k04*k11*k06*k10+k01.*A*k04*k11*k07*k10+k01.*A*k11*k07*k10*k05+k01.*A*k04*k06*k08.*P*k10+k04*k06*k08.*P*k12.*Q*k10;
E3 = k03.*B*k06*k01.*A*k11*k10+k03.*B*k01.*A*k11*k07*k10+k06*k08.*P*k12.*Q*k10*k02+k03.*B*k06*k01.*A*k08.*P*k10+k03.*B*k06*k08.*P*k12.*Q*k10;
E4 = k05*k03.*B*k01.*A*k11*k10+k08.*P*k12.*Q*k10*k02*k04+k05*k08.*P*k12.*Q*k10*k02+k05*k08.*P*k03.*B*k10*k01.*A+k05*k08.*P*k03.*B*k12.*Q*k10;
E5 = k12.*Q*k10*k02*k04*k06+k12.*Q*k07*k10*k02*k04+k12.*Q*k07*k10*k02*k05+k07*k10*k05*k03.*B*k01.*A+k12.*Q*k07*k10*k05*k03.*B;
E6 = k09.*I*k12.*Q*k02*k04*k06+k09.*I*k12.*Q*k07*k02*k04+k09.*I*k12.*Q*k07*k02*k05+k09.*I*k07*k05*k03.*B*k01.*A+k09.*I*k12.*Q*k07*k05*k03.*B;
% Denominator terms
D = E1+E2+E3+E4+E5+E6;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
E6 = E6./D;
% Reaction rate 
v = +k07.*E4-k08.*P.*E5;