function [v,E1,E2] = r_R_GTHOr1(X,K)
% Metabolites definition 
A = X(1);
B = X(2);
P = X(3);
Q = X(4);
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
k13 = K(13);
k14 = K(14);
%  Numerator terms
E1 = k02*k13*k04*k06.*P*k08*k10+k02*k13*k04*k11*k06.*P*k08+k02*k13*k04*k11*k06.*P*k09+k02*k13*k04*k11*k09*k07.*B+k02*k13*k11*k09*k07.*B*k05+k02*k04*k06.*P*k08*k10*k12.*Q+k13*k11*k09*k07.*B*k05*k03;
E2 = k01.*A*k04*k13*k06.*P*k08*k10+k01.*A*k04*k13*k06.*P*k11*k08+k01.*A*k04*k13*k06.*P*k11*k09+k01.*A*k04*k13*k11*k09*k07.*B+k01.*A*k13*k11*k09*k07.*B*k05+k01.*A*k04*k06.*P*k08*k10*k12.*Q+k04*k06.*P*k08*k10*k12.*Q*k14.*Q;
E3 = k03*k06.*P*k01.*A*k08*k13*k10+k03*k06.*P*k01.*A*k08*k13*k11+k03*k06.*P*k01.*A*k13*k11*k09+k03*k01.*A*k13*k11*k09*k07.*B+k06.*P*k08*k10*k12.*Q*k14.*Q*k02+k03*k06.*P*k01.*A*k08*k10*k12.*Q+k03*k06.*P*k08*k10*k12.*Q*k14.*Q;
E4 = k05*k08*k03*k10*k01.*A*k13+k05*k08*k03*k01.*A*k13*k11+k05*k03*k01.*A*k13*k11*k09+k08*k10*k12.*Q*k14.*Q*k02*k04+k05*k08*k10*k12.*Q*k14.*Q*k02+k05*k08*k03*k10*k01.*A*k12.*Q+k05*k08*k03*k10*k12.*Q*k14.*Q;
E5 = k07.*B*k10*k05*k03*k01.*A*k13+k07.*B*k05*k03*k01.*A*k13*k11+k10*k12.*Q*k14.*Q*k02*k04*k06.*P+k07.*B*k10*k12.*Q*k14.*Q*k02*k04+k07.*B*k10*k05*k12.*Q*k14.*Q*k02+k07.*B*k10*k05*k12.*Q*k03*k01.*A+k07.*B*k10*k05*k12.*Q*k03*k14.*Q;
E6 = k09*k07.*B*k05*k03*k01.*A*k13+k12.*Q*k14.*Q*k02*k04*k06.*P*k08+k09*k12.*Q*k14.*Q*k02*k04*k06.*P+k09*k12.*Q*k07.*B*k14.*Q*k02*k04+k09*k12.*Q*k07.*B*k14.*Q*k05*k02+k09*k12.*Q*k07.*B*k05*k03*k01.*A+k09*k12.*Q*k07.*B*k14.*Q*k05*k03;
E7 = k14.*Q*k02*k04*k06.*P*k08*k10+k14.*Q*k11*k02*k04*k06.*P*k08+k14.*Q*k11*k02*k09*k04*k06.*P+k14.*Q*k11*k02*k09*k04*k07.*B+k14.*Q*k11*k02*k09*k07.*B*k05+k11*k09*k07.*B*k05*k03*k01.*A+k14.*Q*k11*k09*k07.*B*k05*k03;
% Denominator terms
D = E1+E2+E3+E4+E5+E6+E7;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
E6 = E6./D;
E7 = E7./D;
% Reaction rate 
v = +k05.*E3-k06.*P.*E4;