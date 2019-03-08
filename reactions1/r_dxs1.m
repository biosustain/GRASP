function [v,E1,E2] = r_dxs1(X,K)
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
k11 = K(11);
k12 = K(12);
%  Numerator terms
E1 = k02*k04*k11*k06+k02*k04*k11*k08+k02*k04*k11*k09+k02*k04*k06*k10.*P+k02*k04*k08*k10.*P+k02*k11*k06*k07.*A+k02*k11*k09*k07.*A+k02*k06*k07.*A*k10.*P+k04*k11*k08*k05.*B+k04*k11*k09*k05.*B+k04*k08*k05.*B*k10.*P+k11*k09*k05.*B*k07.*A;
E2 = k01.*A*k06*k04*k11+k01.*A*k04*k11*k08+k01.*A*k04*k11*k09+k01.*A*k06*k04*k10.*P+k01.*A*k04*k08*k10.*P+k01.*A*k06*k11*k07.*A+k01.*A*k11*k09*k07.*A+k01.*A*k06*k07.*A*k10.*P+k06*k07.*A*k03.*B*k11+k06*k10.*P*k12.*Q*k04+k06*k07.*A*k10.*P*k03.*B+k06*k07.*A*k10.*P*k12.*Q;
E3 = k03.*B*k02*k11*k06+k03.*B*k08*k02*k11+k03.*B*k02*k11*k09+k03.*B*k02*k06*k10.*P+k03.*B*k08*k02*k10.*P+k08*k05.*B*k01.*A*k11+k08*k10.*P*k12.*Q*k02+k08*k05.*B*k10.*P*k01.*A+k03.*B*k08*k11*k05.*B+k03.*B*k11*k09*k05.*B+k03.*B*k08*k05.*B*k10.*P+k08*k05.*B*k10.*P*k12.*Q;
E4 = k05.*B*k01.*A*k04*k11+k07.*A*k03.*B*k02*k11+k10.*P*k12.*Q*k02*k04+k05.*B*k10.*P*k01.*A*k04+k07.*A*k10.*P*k03.*B*k02+k05.*B*k07.*A*k01.*A*k11+k07.*A*k10.*P*k12.*Q*k02+k05.*B*k07.*A*k10.*P*k01.*A+k05.*B*k07.*A*k03.*B*k11+k05.*B*k10.*P*k12.*Q*k04+k05.*B*k07.*A*k10.*P*k03.*B+k05.*B*k07.*A*k10.*P*k12.*Q;
E5 = k12.*Q*k02*k04*k06+k12.*Q*k02*k04*k08+k12.*Q*k09*k02*k04+k09*k05.*B*k01.*A*k04+k09*k07.*A*k03.*B*k02+k12.*Q*k02*k06*k07.*A+k12.*Q*k09*k02*k07.*A+k09*k05.*B*k07.*A*k01.*A+k12.*Q*k04*k08*k05.*B+k12.*Q*k09*k04*k05.*B+k09*k05.*B*k07.*A*k03.*B+k12.*Q*k09*k05.*B*k07.*A;
% Denominator terms
D = E1+E2+E3+E4+E5;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
% Reaction rate 
v = +k09.*E4-k10.*P.*E5;