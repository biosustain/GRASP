function [v,E1,E2] = r_PGK1(X,K)
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
k13 = K(13);
k14 = K(14);
k15 = K(15);
k16 = K(16);
k17 = K(17);
k18 = K(18);
%  Numerator terms
E1 = k02*k04*k15*k17*k06*k10+k02*k04*k15*k17*k06*k11+k02*k04*k15*k17*k06*k13+k02*k04*k15*k17*k08*k10+k02*k04*k15*k17*k08*k11+k02*k04*k15*k17*k08*k13+k02*k04*k15*k17*k11*k09+k02*k04*k15*k17*k13*k09+k02*k04*k15*k06*k10*k14.*Q+k02*k04*k15*k06*k11*k14.*Q+k02*k04*k15*k08*k10*k14.*Q+k02*k04*k15*k08*k11*k14.*Q+k02*k04*k15*k11*k09*k14.*Q+k02*k04*k17*k06*k10*k12.*P+k02*k04*k17*k06*k13*k12.*P+k02*k04*k17*k08*k10*k12.*P+k02*k04*k17*k08*k13*k12.*P+k02*k04*k17*k13*k09*k12.*P+k02*k04*k06*k10*k12.*P*k14.*Q+k02*k04*k08*k10*k12.*P*k14.*Q+k02*k15*k17*k06*k07.*A*k10+k02*k15*k17*k06*k11*k07.*A+k02*k15*k17*k06*k13*k07.*A+k02*k15*k17*k11*k09*k07.*A+k02*k15*k17*k13*k09*k07.*A+k02*k15*k06*k07.*A*k10*k14.*Q+k02*k15*k06*k11*k07.*A*k14.*Q+k02*k15*k11*k09*k14.*Q*k07.*A+k02*k17*k06*k07.*A*k10*k12.*P+k02*k17*k06*k13*k07.*A*k12.*P+k02*k17*k13*k09*k12.*P*k07.*A+k02*k06*k07.*A*k10*k12.*P*k14.*Q+k04*k15*k17*k08*k05.*B*k10+k04*k15*k17*k08*k11*k05.*B+k04*k15*k17*k08*k13*k05.*B+k04*k15*k17*k11*k09*k05.*B+k04*k15*k17*k13*k09*k05.*B+k04*k15*k08*k05.*B*k10*k14.*Q+k04*k15*k08*k11*k05.*B*k14.*Q+k04*k15*k11*k09*k14.*Q*k05.*B+k04*k17*k08*k05.*B*k10*k12.*P+k04*k17*k08*k13*k05.*B*k12.*P+k04*k17*k13*k09*k12.*P*k05.*B+k04*k08*k05.*B*k10*k12.*P*k14.*Q+k15*k17*k11*k09*k05.*B*k07.*A+k15*k17*k13*k09*k05.*B*k07.*A+k15*k11*k09*k14.*Q*k05.*B*k07.*A+k17*k13*k09*k12.*P*k05.*B*k07.*A;
E2 = k01.*A*k06*k04*k15*k17*k10+k01.*A*k06*k04*k15*k17*k11+k01.*A*k06*k04*k15*k17*k13+k01.*A*k04*k15*k17*k08*k10+k01.*A*k04*k15*k17*k08*k11+k01.*A*k04*k15*k17*k08*k13+k01.*A*k04*k15*k17*k11*k09+k01.*A*k04*k15*k17*k13*k09+k01.*A*k06*k04*k15*k10*k14.*Q+k01.*A*k06*k04*k15*k11*k14.*Q+k01.*A*k04*k15*k08*k10*k14.*Q+k01.*A*k04*k15*k08*k11*k14.*Q+k01.*A*k04*k15*k11*k09*k14.*Q+k01.*A*k06*k04*k17*k10*k12.*P+k01.*A*k06*k04*k17*k13*k12.*P+k01.*A*k04*k17*k08*k10*k12.*P+k01.*A*k04*k17*k08*k13*k12.*P+k01.*A*k04*k17*k13*k09*k12.*P+k01.*A*k06*k04*k10*k12.*P*k14.*Q+k01.*A*k04*k08*k10*k12.*P*k14.*Q+k01.*A*k06*k15*k17*k07.*A*k10+k01.*A*k06*k15*k17*k07.*A*k11+k01.*A*k06*k15*k17*k07.*A*k13+k01.*A*k15*k17*k11*k09*k07.*A+k01.*A*k15*k17*k13*k09*k07.*A+k01.*A*k06*k15*k07.*A*k10*k14.*Q+k01.*A*k06*k15*k07.*A*k11*k14.*Q+k01.*A*k15*k11*k09*k14.*Q*k07.*A+k01.*A*k06*k17*k07.*A*k10*k12.*P+k01.*A*k06*k17*k07.*A*k13*k12.*P+k01.*A*k17*k13*k09*k12.*P*k07.*A+k01.*A*k06*k07.*A*k10*k12.*P*k14.*Q+k06*k07.*A*k10*k03.*B*k15*k17+k06*k07.*A*k03.*B*k15*k17*k11+k06*k07.*A*k03.*B*k15*k17*k13+k06*k10*k12.*P*k16.*Q*k04*k17+k06*k10*k14.*Q*k18.*P*k04*k15+k06*k07.*A*k10*k03.*B*k14.*Q*k15+k06*k07.*A*k03.*B*k15*k11*k14.*Q+k06*k10*k12.*P*k14.*Q*k16.*Q*k04+k06*k07.*A*k10*k03.*B*k12.*P*k17+k06*k07.*A*k03.*B*k17*k13*k12.*P+k06*k10*k12.*P*k14.*Q*k18.*P*k04+k06*k07.*A*k10*k03.*B*k12.*P*k14.*Q+k06*k07.*A*k10*k12.*P*k16.*Q*k17+k06*k07.*A*k10*k14.*Q*k18.*P*k15+k06*k07.*A*k10*k12.*P*k14.*Q*k16.*Q+k06*k07.*A*k10*k12.*P*k14.*Q*k18.*P;
E3 = k03.*B*k02*k15*k17*k06*k10+k03.*B*k02*k15*k17*k06*k11+k03.*B*k02*k15*k17*k06*k13+k03.*B*k08*k02*k15*k17*k10+k03.*B*k08*k02*k15*k17*k11+k03.*B*k08*k02*k15*k17*k13+k03.*B*k02*k15*k17*k11*k09+k03.*B*k02*k15*k17*k13*k09+k03.*B*k02*k15*k06*k10*k14.*Q+k03.*B*k02*k15*k06*k11*k14.*Q+k03.*B*k08*k02*k15*k10*k14.*Q+k03.*B*k08*k02*k15*k11*k14.*Q+k03.*B*k02*k15*k11*k09*k14.*Q+k03.*B*k02*k17*k06*k10*k12.*P+k03.*B*k02*k17*k06*k13*k12.*P+k03.*B*k08*k02*k17*k10*k12.*P+k03.*B*k08*k02*k17*k13*k12.*P+k03.*B*k02*k17*k13*k09*k12.*P+k03.*B*k02*k06*k10*k12.*P*k14.*Q+k03.*B*k08*k02*k10*k12.*P*k14.*Q+k08*k05.*B*k10*k01.*A*k15*k17+k08*k05.*B*k01.*A*k15*k17*k11+k08*k05.*B*k01.*A*k15*k17*k13+k08*k10*k12.*P*k16.*Q*k02*k17+k08*k10*k14.*Q*k18.*P*k02*k15+k08*k05.*B*k10*k01.*A*k14.*Q*k15+k08*k05.*B*k01.*A*k15*k11*k14.*Q+k08*k10*k12.*P*k14.*Q*k16.*Q*k02+k08*k05.*B*k10*k01.*A*k12.*P*k17+k08*k05.*B*k01.*A*k17*k13*k12.*P+k08*k10*k12.*P*k14.*Q*k18.*P*k02+k08*k05.*B*k10*k01.*A*k12.*P*k14.*Q+k03.*B*k08*k15*k17*k05.*B*k10+k03.*B*k08*k15*k17*k05.*B*k11+k03.*B*k08*k15*k17*k05.*B*k13+k03.*B*k15*k17*k11*k09*k05.*B+k03.*B*k15*k17*k13*k09*k05.*B+k03.*B*k08*k15*k05.*B*k10*k14.*Q+k03.*B*k08*k15*k05.*B*k11*k14.*Q+k03.*B*k15*k11*k09*k14.*Q*k05.*B+k03.*B*k08*k17*k05.*B*k10*k12.*P+k03.*B*k08*k17*k05.*B*k13*k12.*P+k03.*B*k17*k13*k09*k12.*P*k05.*B+k03.*B*k08*k05.*B*k10*k12.*P*k14.*Q+k08*k05.*B*k10*k12.*P*k16.*Q*k17+k08*k05.*B*k10*k14.*Q*k18.*P*k15+k08*k05.*B*k10*k12.*P*k14.*Q*k16.*Q+k08*k05.*B*k10*k12.*P*k14.*Q*k18.*P;
E4 = k05.*B*k10*k01.*A*k04*k15*k17+k05.*B*k01.*A*k04*k15*k17*k11+k05.*B*k01.*A*k04*k15*k17*k13+k07.*A*k10*k03.*B*k02*k15*k17+k07.*A*k03.*B*k02*k15*k17*k11+k07.*A*k03.*B*k02*k15*k17*k13+k10*k12.*P*k16.*Q*k02*k04*k17+k10*k14.*Q*k18.*P*k02*k04*k15+k05.*B*k10*k01.*A*k14.*Q*k04*k15+k05.*B*k01.*A*k04*k15*k11*k14.*Q+k07.*A*k10*k03.*B*k14.*Q*k02*k15+k07.*A*k03.*B*k02*k15*k11*k14.*Q+k10*k12.*P*k14.*Q*k16.*Q*k02*k04+k05.*B*k10*k01.*A*k12.*P*k04*k17+k05.*B*k01.*A*k04*k17*k13*k12.*P+k07.*A*k10*k03.*B*k12.*P*k02*k17+k07.*A*k03.*B*k02*k17*k13*k12.*P+k10*k12.*P*k14.*Q*k18.*P*k02*k04+k05.*B*k10*k01.*A*k12.*P*k14.*Q*k04+k07.*A*k10*k03.*B*k12.*P*k14.*Q*k02+k05.*B*k07.*A*k10*k01.*A*k15*k17+k05.*B*k07.*A*k01.*A*k15*k17*k11+k05.*B*k07.*A*k01.*A*k15*k17*k13+k07.*A*k10*k12.*P*k16.*Q*k02*k17+k07.*A*k10*k14.*Q*k18.*P*k02*k15+k05.*B*k07.*A*k10*k01.*A*k14.*Q*k15+k05.*B*k07.*A*k01.*A*k15*k11*k14.*Q+k07.*A*k10*k12.*P*k14.*Q*k16.*Q*k02+k05.*B*k07.*A*k10*k01.*A*k12.*P*k17+k05.*B*k07.*A*k01.*A*k17*k13*k12.*P+k07.*A*k10*k12.*P*k14.*Q*k18.*P*k02+k05.*B*k07.*A*k10*k01.*A*k12.*P*k14.*Q+k05.*B*k07.*A*k10*k03.*B*k15*k17+k05.*B*k07.*A*k03.*B*k15*k17*k11+k05.*B*k07.*A*k03.*B*k15*k17*k13+k05.*B*k10*k12.*P*k16.*Q*k04*k17+k05.*B*k10*k14.*Q*k18.*P*k04*k15+k05.*B*k07.*A*k10*k03.*B*k14.*Q*k15+k05.*B*k07.*A*k03.*B*k15*k11*k14.*Q+k05.*B*k10*k12.*P*k14.*Q*k16.*Q*k04+k05.*B*k07.*A*k10*k03.*B*k12.*P*k17+k05.*B*k07.*A*k03.*B*k17*k13*k12.*P+k05.*B*k10*k12.*P*k14.*Q*k18.*P*k04+k05.*B*k07.*A*k10*k03.*B*k12.*P*k14.*Q+k05.*B*k07.*A*k10*k12.*P*k16.*Q*k17+k05.*B*k07.*A*k10*k14.*Q*k18.*P*k15+k05.*B*k07.*A*k10*k12.*P*k14.*Q*k16.*Q+k05.*B*k07.*A*k10*k12.*P*k14.*Q*k18.*P;
E5 = k09*k05.*B*k01.*A*k04*k15*k17+k12.*P*k16.*Q*k02*k04*k17*k06+k14.*Q*k18.*P*k02*k04*k15*k06+k09*k07.*A*k03.*B*k02*k15*k17+k12.*P*k16.*Q*k02*k04*k17*k08+k14.*Q*k18.*P*k02*k04*k15*k08+k09*k12.*P*k16.*Q*k02*k04*k17+k09*k14.*Q*k18.*P*k02*k04*k15+k09*k14.*Q*k05.*B*k01.*A*k04*k15+k12.*P*k14.*Q*k16.*Q*k02*k04*k06+k09*k14.*Q*k07.*A*k03.*B*k02*k15+k12.*P*k14.*Q*k16.*Q*k02*k04*k08+k09*k12.*P*k14.*Q*k16.*Q*k02*k04+k09*k12.*P*k05.*B*k01.*A*k04*k17+k12.*P*k14.*Q*k18.*P*k02*k04*k06+k09*k12.*P*k07.*A*k03.*B*k02*k17+k12.*P*k14.*Q*k18.*P*k02*k04*k08+k09*k12.*P*k14.*Q*k18.*P*k02*k04+k09*k12.*P*k14.*Q*k05.*B*k01.*A*k04+k09*k12.*P*k14.*Q*k07.*A*k03.*B*k02+k09*k05.*B*k07.*A*k01.*A*k15*k17+k12.*P*k16.*Q*k02*k17*k06*k07.*A+k14.*Q*k18.*P*k02*k15*k06*k07.*A+k09*k12.*P*k07.*A*k16.*Q*k02*k17+k09*k14.*Q*k07.*A*k18.*P*k02*k15+k09*k14.*Q*k05.*B*k07.*A*k01.*A*k15+k12.*P*k14.*Q*k16.*Q*k02*k06*k07.*A+k09*k12.*P*k14.*Q*k07.*A*k16.*Q*k02+k09*k12.*P*k05.*B*k07.*A*k01.*A*k17+k12.*P*k14.*Q*k18.*P*k02*k06*k07.*A+k09*k12.*P*k14.*Q*k07.*A*k18.*P*k02+k09*k12.*P*k14.*Q*k05.*B*k07.*A*k01.*A+k09*k05.*B*k07.*A*k03.*B*k15*k17+k12.*P*k16.*Q*k04*k17*k08*k05.*B+k14.*Q*k18.*P*k04*k15*k08*k05.*B+k09*k12.*P*k05.*B*k16.*Q*k04*k17+k09*k14.*Q*k05.*B*k18.*P*k04*k15+k09*k14.*Q*k05.*B*k07.*A*k03.*B*k15+k12.*P*k14.*Q*k16.*Q*k04*k08*k05.*B+k09*k12.*P*k14.*Q*k05.*B*k16.*Q*k04+k09*k12.*P*k05.*B*k07.*A*k03.*B*k17+k12.*P*k14.*Q*k18.*P*k04*k08*k05.*B+k09*k12.*P*k14.*Q*k05.*B*k18.*P*k04+k09*k12.*P*k14.*Q*k05.*B*k07.*A*k03.*B+k09*k12.*P*k05.*B*k07.*A*k16.*Q*k17+k09*k14.*Q*k05.*B*k07.*A*k18.*P*k15+k09*k12.*P*k14.*Q*k05.*B*k07.*A*k16.*Q+k09*k12.*P*k14.*Q*k05.*B*k07.*A*k18.*P;
E6 = k16.*Q*k02*k04*k17*k06*k10+k16.*Q*k11*k02*k04*k17*k06+k16.*Q*k02*k04*k17*k06*k13+k16.*Q*k02*k04*k17*k08*k10+k16.*Q*k11*k02*k04*k17*k08+k16.*Q*k02*k04*k17*k08*k13+k16.*Q*k11*k02*k04*k17*k09+k16.*Q*k02*k04*k17*k13*k09+k16.*Q*k02*k04*k06*k10*k14.*Q+k16.*Q*k11*k02*k04*k14.*Q*k06+k16.*Q*k02*k04*k08*k10*k14.*Q+k16.*Q*k11*k02*k04*k14.*Q*k08+k16.*Q*k11*k02*k04*k09*k14.*Q+k11*k09*k05.*B*k01.*A*k04*k17+k11*k14.*Q*k18.*P*k02*k04*k06+k11*k09*k07.*A*k03.*B*k02*k17+k11*k14.*Q*k18.*P*k02*k04*k08+k11*k09*k14.*Q*k18.*P*k02*k04+k11*k09*k14.*Q*k05.*B*k01.*A*k04+k11*k09*k14.*Q*k07.*A*k03.*B*k02+k16.*Q*k02*k17*k06*k07.*A*k10+k16.*Q*k11*k02*k17*k06*k07.*A+k16.*Q*k02*k17*k06*k13*k07.*A+k16.*Q*k11*k02*k17*k09*k07.*A+k16.*Q*k02*k17*k13*k09*k07.*A+k16.*Q*k02*k06*k07.*A*k10*k14.*Q+k16.*Q*k11*k02*k14.*Q*k06*k07.*A+k16.*Q*k11*k02*k09*k14.*Q*k07.*A+k11*k09*k05.*B*k07.*A*k01.*A*k17+k11*k14.*Q*k18.*P*k02*k06*k07.*A+k11*k09*k14.*Q*k07.*A*k18.*P*k02+k11*k09*k14.*Q*k05.*B*k07.*A*k01.*A+k16.*Q*k04*k17*k08*k05.*B*k10+k16.*Q*k11*k04*k17*k08*k05.*B+k16.*Q*k04*k17*k08*k13*k05.*B+k16.*Q*k11*k04*k17*k09*k05.*B+k16.*Q*k04*k17*k13*k09*k05.*B+k16.*Q*k04*k08*k05.*B*k10*k14.*Q+k16.*Q*k11*k04*k14.*Q*k08*k05.*B+k16.*Q*k11*k04*k09*k14.*Q*k05.*B+k11*k09*k05.*B*k07.*A*k03.*B*k17+k11*k14.*Q*k18.*P*k04*k08*k05.*B+k11*k09*k14.*Q*k05.*B*k18.*P*k04+k11*k09*k14.*Q*k05.*B*k07.*A*k03.*B+k16.*Q*k11*k17*k09*k05.*B*k07.*A+k16.*Q*k17*k13*k09*k05.*B*k07.*A+k16.*Q*k11*k09*k14.*Q*k05.*B*k07.*A+k11*k09*k14.*Q*k05.*B*k07.*A*k18.*P;
E7 = k18.*P*k02*k04*k15*k06*k10+k18.*P*k02*k04*k15*k06*k11+k18.*P*k13*k02*k04*k15*k06+k18.*P*k02*k04*k15*k08*k10+k18.*P*k02*k04*k15*k08*k11+k18.*P*k13*k02*k04*k15*k08+k18.*P*k02*k04*k15*k11*k09+k18.*P*k13*k02*k04*k15*k09+k13*k09*k05.*B*k01.*A*k04*k15+k13*k12.*P*k16.*Q*k02*k04*k06+k13*k09*k07.*A*k03.*B*k02*k15+k13*k12.*P*k16.*Q*k02*k04*k08+k13*k09*k12.*P*k16.*Q*k02*k04+k18.*P*k02*k04*k06*k10*k12.*P+k18.*P*k13*k02*k04*k12.*P*k06+k18.*P*k02*k04*k08*k10*k12.*P+k18.*P*k13*k02*k04*k12.*P*k08+k18.*P*k13*k02*k04*k09*k12.*P+k13*k09*k12.*P*k05.*B*k01.*A*k04+k13*k09*k12.*P*k07.*A*k03.*B*k02+k18.*P*k02*k15*k06*k07.*A*k10+k18.*P*k02*k15*k06*k11*k07.*A+k18.*P*k13*k02*k15*k06*k07.*A+k18.*P*k02*k15*k11*k09*k07.*A+k18.*P*k13*k02*k15*k09*k07.*A+k13*k09*k05.*B*k07.*A*k01.*A*k15+k13*k12.*P*k16.*Q*k02*k06*k07.*A+k13*k09*k12.*P*k07.*A*k16.*Q*k02+k18.*P*k02*k06*k07.*A*k10*k12.*P+k18.*P*k13*k02*k12.*P*k06*k07.*A+k18.*P*k13*k02*k09*k12.*P*k07.*A+k13*k09*k12.*P*k05.*B*k07.*A*k01.*A+k18.*P*k04*k15*k08*k05.*B*k10+k18.*P*k04*k15*k08*k11*k05.*B+k18.*P*k13*k04*k15*k08*k05.*B+k18.*P*k04*k15*k11*k09*k05.*B+k18.*P*k13*k04*k15*k09*k05.*B+k13*k09*k05.*B*k07.*A*k03.*B*k15+k13*k12.*P*k16.*Q*k04*k08*k05.*B+k13*k09*k12.*P*k05.*B*k16.*Q*k04+k18.*P*k04*k08*k05.*B*k10*k12.*P+k18.*P*k13*k04*k12.*P*k08*k05.*B+k18.*P*k13*k04*k09*k12.*P*k05.*B+k13*k09*k12.*P*k05.*B*k07.*A*k03.*B+k18.*P*k15*k11*k09*k05.*B*k07.*A+k18.*P*k13*k15*k09*k05.*B*k07.*A+k13*k09*k12.*P*k05.*B*k07.*A*k16.*Q+k18.*P*k13*k09*k12.*P*k05.*B*k07.*A;
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
v = +k17.*E7-k18.*P.*E1+k11.*E5-k12.*P.*E6;