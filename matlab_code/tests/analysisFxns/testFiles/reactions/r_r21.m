function [v,E1,E2] = r_r21(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
P1 = X(3,:);
P2 = X(4,:);
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
E1 = k02*k05*k08*k11+k02*k05*k08*k10+k02*k05*k11*k09+k02*k08*k11*k04+k02*k08*k04*k10+k02*k11*k04*k09+k05*k08*k11*k03+k05*k08*k03*k10+k05*k11*k03*k09;
E2 = k01.*A*k05*k08*k11+k01.*A*k05*k08*k10+k01.*A*k05*k11*k09+k01.*A*k04*k08*k11+k01.*A*k04*k08*k10+k01.*A*k04*k11*k09+k04*k06.*P1*k08*k11+k04*k06.*P1*k08*k10+k04*k06.*P1*k11*k09;
E3 = k06.*P1*k02*k08*k11+k06.*P1*k02*k08*k10+k06.*P1*k02*k11*k09+k03*k01.*A*k08*k11+k03*k01.*A*k08*k10+k03*k01.*A*k11*k09+k06.*P1*k03*k08*k11+k06.*P1*k03*k08*k10+k06.*P1*k03*k11*k09;
E4 = k07.*B*k02*k05*k11+k07.*B*k10*k02*k05+k10*k12.*P2*k02*k05+k07.*B*k02*k11*k04+k07.*B*k10*k02*k04+k10*k12.*P2*k02*k04+k07.*B*k05*k11*k03+k07.*B*k10*k05*k03+k10*k12.*P2*k05*k03;
E5 = k12.*P2*k02*k05*k08+k09*k07.*B*k02*k05+k12.*P2*k09*k02*k05+k12.*P2*k02*k08*k04+k09*k07.*B*k02*k04+k12.*P2*k09*k02*k04+k12.*P2*k05*k08*k03+k09*k07.*B*k05*k03+k12.*P2*k09*k05*k03;
% Denominator terms
D = E1+E2+E3+E4+E5;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
% Reaction rate 
v = +k05.*E3-k06.*P1.*E1;