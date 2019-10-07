function [v,E1,E2] = r_r11(X,K)
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
E1 = k02*k04*k11*k06*k08+k02*k04*k11*k06*k09+k02*k04*k11*k09*k07+k02*k04*k06*k08*k10.*P+k04*k11*k09*k07*k05.*B;
E2 = k01.*A*k06*k04*k11*k08+k01.*A*k06*k04*k11*k09+k01.*A*k04*k11*k09*k07+k01.*A*k06*k04*k08*k10.*P+k06*k08*k10.*P*k12.*Q*k04;
E3 = k03.*I*k02*k11*k06*k08+k03.*I*k02*k11*k06*k09+k03.*I*k02*k11*k09*k07+k03.*I*k02*k06*k08*k10.*P+k03.*I*k11*k09*k07*k05.*B;
E4 = k05.*B*k08*k01.*A*k04*k11+k05.*B*k01.*A*k04*k11*k09+k08*k10.*P*k12.*Q*k02*k04+k05.*B*k08*k01.*A*k10.*P*k04+k05.*B*k08*k10.*P*k12.*Q*k04;
E5 = k07*k05.*B*k01.*A*k04*k11+k10.*P*k12.*Q*k02*k04*k06+k07*k10.*P*k12.*Q*k02*k04+k07*k10.*P*k05.*B*k01.*A*k04+k07*k10.*P*k05.*B*k12.*Q*k04;
E6 = k12.*Q*k02*k04*k06*k08+k12.*Q*k09*k02*k04*k06+k12.*Q*k09*k02*k04*k07+k09*k07*k05.*B*k01.*A*k04+k12.*Q*k09*k04*k07*k05.*B;
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
v = +k09.*E5-k10.*P.*E6;