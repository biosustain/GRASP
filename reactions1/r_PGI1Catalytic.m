function [v,E1,E2] = r_PGI1Catalytic(X,K)
% Metabolites definition 
A = X(1,:);
I1 = X(2,:);
P = X(3,:);
Z1 = X(4,:);
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
E1 = k02*k04*k11*k07+k02*k04*k11*k06+k02*k04*k07*k10+k02*k04*k10*k06+k02*k11*k07*k05.*A+k02*k07*k10*k05.*A+k04*k11*k09*k07+k04*k11*k06*k08.*I1+k04*k11*k06*k09+k04*k06*k08.*I1*k10+k11*k09*k07*k05.*A;
E2 = k01.*A*k07*k04*k11+k01.*A*k04*k11*k06+k01.*A*k07*k10*k04+k01.*A*k10*k04*k06+k01.*A*k07*k11*k05.*A+k01.*A*k07*k10*k05.*A+k07*k10*k12.*P*k04+k07*k05.*A*k03.*Z1*k11+k10*k12.*P*k04*k06+k07*k10*k05.*A*k03.*Z1+k07*k10*k05.*A*k12.*P;
E3 = k03.*Z1*k02*k11*k07+k03.*Z1*k06*k02*k11+k03.*Z1*k02*k07*k10+k03.*Z1*k06*k02*k10+k06*k08.*I1*k01.*A*k11+k06*k08.*I1*k01.*A*k10+k03.*Z1*k11*k09*k07+k03.*Z1*k06*k11*k08.*I1+k03.*Z1*k06*k11*k09+k03.*Z1*k06*k08.*I1*k10+k06*k08.*I1*k10*k12.*P;
E4 = k08.*I1*k01.*A*k04*k11+k05.*A*k03.*Z1*k02*k11+k08.*I1*k01.*A*k10*k04+k05.*A*k03.*Z1*k02*k10+k08.*I1*k05.*A*k01.*A*k11+k08.*I1*k05.*A*k01.*A*k10+k08.*I1*k10*k12.*P*k04+k08.*I1*k05.*A*k03.*Z1*k11+k05.*A*k03.*Z1*k11*k09+k08.*I1*k05.*A*k10*k03.*Z1+k08.*I1*k05.*A*k10*k12.*P;
E5 = k12.*P*k02*k04*k07+k12.*P*k02*k04*k06+k09*k01.*A*k07*k04+k09*k01.*A*k04*k06+k12.*P*k02*k07*k05.*A+k09*k01.*A*k07*k05.*A+k12.*P*k09*k04*k07+k12.*P*k04*k06*k08.*I1+k12.*P*k09*k04*k06+k09*k07*k05.*A*k03.*Z1+k12.*P*k09*k07*k05.*A;
% Denominator terms
D = E1+E2+E3+E4+E5;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
% Reaction rate 
v = +k11.*E5-k12.*P.*E1;