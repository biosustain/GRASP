function [v,E1,E2] = r_ASMT1(X,K)
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
k13 = K(13);
k14 = K(14);
k15 = K(15);
k16 = K(16);
k17 = K(17);
k18 = K(18);
k19 = K(19);
k20 = K(20);
k21 = K(21);
k22 = K(22);
%  Numerator terms
E1 = k02*k04*k06*k19*k21*k10*k08*k14+k02*k04*k06*k19*k21*k10*k08*k15+k02*k04*k06*k19*k21*k10*k08*k17+k02*k04*k06*k19*k21*k08*k12*k14+k02*k04*k06*k19*k21*k08*k12*k15+k02*k04*k06*k19*k21*k08*k12*k17+k02*k04*k06*k19*k21*k08*k15*k13+k02*k04*k06*k19*k21*k08*k17*k13+k02*k04*k06*k19*k10*k08*k14*k18.*Q+k02*k04*k06*k19*k10*k08*k15*k18.*Q+k02*k04*k06*k19*k08*k12*k14*k18.*Q+k02*k04*k06*k19*k08*k12*k15*k18.*Q+k02*k04*k06*k19*k08*k15*k13*k18.*Q+k02*k04*k06*k21*k10*k08*k14*k16.*P+k02*k04*k06*k21*k10*k08*k17*k16.*P+k02*k04*k06*k21*k08*k12*k14*k16.*P+k02*k04*k06*k21*k08*k12*k17*k16.*P+k02*k04*k06*k21*k08*k17*k13*k16.*P+k02*k04*k06*k10*k08*k14*k16.*P*k18.*Q+k02*k04*k06*k08*k12*k14*k16.*P*k18.*Q+k02*k04*k19*k21*k10*k11.*A*k14*k08+k02*k04*k19*k21*k10*k15*k11.*A*k08+k02*k04*k19*k21*k10*k17*k11.*A*k08+k02*k04*k19*k21*k15*k13*k11.*A*k08+k02*k04*k19*k21*k17*k13*k11.*A*k08+k02*k04*k19*k10*k11.*A*k14*k08*k18.*Q+k02*k04*k19*k10*k15*k11.*A*k18.*Q*k08+k02*k04*k19*k15*k13*k18.*Q*k11.*A*k08+k02*k04*k21*k10*k11.*A*k14*k08*k16.*P+k02*k04*k21*k10*k17*k11.*A*k16.*P*k08+k02*k04*k21*k17*k13*k16.*P*k11.*A*k08+k02*k04*k10*k11.*A*k14*k08*k16.*P*k18.*Q+k04*k06*k19*k21*k08*k12*k09.*B*k14+k04*k06*k19*k21*k08*k12*k15*k09.*B+k04*k06*k19*k21*k08*k12*k17*k09.*B+k04*k06*k19*k21*k08*k15*k13*k09.*B+k04*k06*k19*k21*k08*k17*k13*k09.*B+k04*k06*k19*k08*k12*k09.*B*k14*k18.*Q+k04*k06*k19*k08*k12*k15*k09.*B*k18.*Q+k04*k06*k19*k08*k15*k13*k18.*Q*k09.*B+k04*k06*k21*k08*k12*k09.*B*k14*k16.*P+k04*k06*k21*k08*k12*k17*k09.*B*k16.*P+k04*k06*k21*k08*k17*k13*k16.*P*k09.*B+k04*k06*k08*k12*k09.*B*k14*k16.*P*k18.*Q+k04*k19*k21*k15*k13*k09.*B*k11.*A*k08+k04*k19*k21*k17*k13*k09.*B*k11.*A*k08+k04*k19*k15*k13*k18.*Q*k09.*B*k11.*A*k08+k04*k21*k17*k13*k16.*P*k09.*B*k11.*A*k08;
E2 = k01.*A*k10*k04*k06*k19*k21*k14*k08+k01.*A*k10*k04*k06*k19*k21*k08*k15+k01.*A*k10*k04*k06*k19*k21*k08*k17+k01.*A*k04*k06*k19*k21*k08*k12*k14+k01.*A*k04*k06*k19*k21*k08*k12*k15+k01.*A*k04*k06*k19*k21*k08*k12*k17+k01.*A*k04*k06*k19*k21*k08*k15*k13+k01.*A*k04*k06*k19*k21*k08*k17*k13+k01.*A*k10*k04*k06*k19*k14*k08*k18.*Q+k01.*A*k10*k04*k06*k19*k08*k15*k18.*Q+k01.*A*k04*k06*k19*k08*k12*k14*k18.*Q+k01.*A*k04*k06*k19*k08*k12*k15*k18.*Q+k01.*A*k04*k06*k19*k08*k15*k13*k18.*Q+k01.*A*k10*k04*k06*k21*k14*k08*k16.*P+k01.*A*k10*k04*k06*k21*k08*k17*k16.*P+k01.*A*k04*k06*k21*k08*k12*k14*k16.*P+k01.*A*k04*k06*k21*k08*k12*k17*k16.*P+k01.*A*k04*k06*k21*k08*k17*k13*k16.*P+k01.*A*k10*k04*k06*k14*k08*k16.*P*k18.*Q+k01.*A*k04*k06*k08*k12*k14*k16.*P*k18.*Q+k01.*A*k10*k04*k19*k21*k11.*A*k14*k08+k01.*A*k10*k04*k19*k21*k11.*A*k15*k08+k01.*A*k10*k04*k19*k21*k11.*A*k17*k08+k01.*A*k04*k19*k21*k15*k13*k11.*A*k08+k01.*A*k04*k19*k21*k17*k13*k11.*A*k08+k01.*A*k10*k04*k19*k11.*A*k14*k08*k18.*Q+k01.*A*k10*k04*k19*k11.*A*k15*k08*k18.*Q+k01.*A*k04*k19*k15*k13*k18.*Q*k11.*A*k08+k01.*A*k10*k04*k21*k11.*A*k14*k08*k16.*P+k01.*A*k10*k04*k21*k11.*A*k17*k08*k16.*P+k01.*A*k04*k21*k17*k13*k16.*P*k11.*A*k08+k01.*A*k10*k04*k11.*A*k14*k08*k16.*P*k18.*Q+k10*k11.*A*k14*k05.*B*k08*k04*k19*k21+k10*k11.*A*k05.*B*k08*k04*k19*k21*k15+k10*k11.*A*k05.*B*k08*k04*k19*k21*k17+k10*k14*k16.*P*k20.*Q*k04*k06*k21*k08+k10*k14*k18.*Q*k22.*P*k04*k06*k19*k08+k10*k11.*A*k14*k05.*B*k08*k18.*Q*k04*k19+k10*k11.*A*k05.*B*k08*k04*k19*k15*k18.*Q+k10*k14*k16.*P*k18.*Q*k20.*Q*k04*k06*k08+k10*k11.*A*k14*k05.*B*k08*k16.*P*k04*k21+k10*k11.*A*k05.*B*k08*k04*k21*k17*k16.*P+k10*k14*k16.*P*k18.*Q*k22.*P*k04*k06*k08+k10*k11.*A*k14*k05.*B*k08*k16.*P*k18.*Q*k04+k10*k11.*A*k14*k08*k16.*P*k20.*Q*k04*k21+k10*k11.*A*k14*k08*k18.*Q*k22.*P*k04*k19+k10*k11.*A*k14*k08*k16.*P*k18.*Q*k20.*Q*k04+k10*k11.*A*k14*k08*k16.*P*k18.*Q*k22.*P*k04;
E3 = k03.*I*k02*k06*k19*k21*k10*k08*k14+k03.*I*k02*k06*k19*k21*k10*k08*k15+k03.*I*k02*k06*k19*k21*k10*k08*k17+k03.*I*k02*k06*k19*k21*k08*k12*k14+k03.*I*k02*k06*k19*k21*k08*k12*k15+k03.*I*k02*k06*k19*k21*k08*k12*k17+k03.*I*k02*k06*k19*k21*k08*k15*k13+k03.*I*k02*k06*k19*k21*k08*k17*k13+k03.*I*k02*k06*k19*k10*k08*k14*k18.*Q+k03.*I*k02*k06*k19*k10*k08*k15*k18.*Q+k03.*I*k02*k06*k19*k08*k12*k14*k18.*Q+k03.*I*k02*k06*k19*k08*k12*k15*k18.*Q+k03.*I*k02*k06*k19*k08*k15*k13*k18.*Q+k03.*I*k02*k06*k21*k10*k08*k14*k16.*P+k03.*I*k02*k06*k21*k10*k08*k17*k16.*P+k03.*I*k02*k06*k21*k08*k12*k14*k16.*P+k03.*I*k02*k06*k21*k08*k12*k17*k16.*P+k03.*I*k02*k06*k21*k08*k17*k13*k16.*P+k03.*I*k02*k06*k10*k08*k14*k16.*P*k18.*Q+k03.*I*k02*k06*k08*k12*k14*k16.*P*k18.*Q+k03.*I*k02*k19*k21*k10*k11.*A*k14*k08+k03.*I*k02*k19*k21*k10*k15*k11.*A*k08+k03.*I*k02*k19*k21*k10*k17*k11.*A*k08+k03.*I*k02*k19*k21*k15*k13*k11.*A*k08+k03.*I*k02*k19*k21*k17*k13*k11.*A*k08+k03.*I*k02*k19*k10*k11.*A*k14*k08*k18.*Q+k03.*I*k02*k19*k10*k15*k11.*A*k18.*Q*k08+k03.*I*k02*k19*k15*k13*k18.*Q*k11.*A*k08+k03.*I*k02*k21*k10*k11.*A*k14*k08*k16.*P+k03.*I*k02*k21*k10*k17*k11.*A*k16.*P*k08+k03.*I*k02*k21*k17*k13*k16.*P*k11.*A*k08+k03.*I*k02*k10*k11.*A*k14*k08*k16.*P*k18.*Q+k03.*I*k06*k19*k21*k08*k12*k09.*B*k14+k03.*I*k06*k19*k21*k08*k12*k15*k09.*B+k03.*I*k06*k19*k21*k08*k12*k17*k09.*B+k03.*I*k06*k19*k21*k08*k15*k13*k09.*B+k03.*I*k06*k19*k21*k08*k17*k13*k09.*B+k03.*I*k06*k19*k08*k12*k09.*B*k14*k18.*Q+k03.*I*k06*k19*k08*k12*k15*k09.*B*k18.*Q+k03.*I*k06*k19*k08*k15*k13*k18.*Q*k09.*B+k03.*I*k06*k21*k08*k12*k09.*B*k14*k16.*P+k03.*I*k06*k21*k08*k12*k17*k09.*B*k16.*P+k03.*I*k06*k21*k08*k17*k13*k16.*P*k09.*B+k03.*I*k06*k08*k12*k09.*B*k14*k16.*P*k18.*Q+k03.*I*k19*k21*k15*k13*k09.*B*k11.*A*k08+k03.*I*k19*k21*k17*k13*k09.*B*k11.*A*k08+k03.*I*k19*k15*k13*k18.*Q*k09.*B*k11.*A*k08+k03.*I*k21*k17*k13*k16.*P*k09.*B*k11.*A*k08;
E4 = k05.*B*k08*k02*k04*k19*k21*k10*k14+k05.*B*k08*k02*k04*k19*k21*k10*k15+k05.*B*k08*k02*k04*k19*k21*k10*k17+k05.*B*k08*k12*k02*k04*k19*k21*k14+k05.*B*k08*k12*k02*k04*k19*k21*k15+k05.*B*k08*k12*k02*k04*k19*k21*k17+k05.*B*k08*k02*k04*k19*k21*k15*k13+k05.*B*k08*k02*k04*k19*k21*k17*k13+k05.*B*k08*k02*k04*k19*k10*k14*k18.*Q+k05.*B*k08*k02*k04*k19*k10*k15*k18.*Q+k05.*B*k08*k12*k02*k04*k19*k14*k18.*Q+k05.*B*k08*k12*k02*k04*k19*k15*k18.*Q+k05.*B*k08*k02*k04*k19*k15*k13*k18.*Q+k05.*B*k08*k02*k04*k21*k10*k14*k16.*P+k05.*B*k08*k02*k04*k21*k10*k17*k16.*P+k05.*B*k08*k12*k02*k04*k21*k14*k16.*P+k05.*B*k08*k12*k02*k04*k21*k17*k16.*P+k05.*B*k08*k02*k04*k21*k17*k13*k16.*P+k05.*B*k08*k02*k04*k10*k14*k16.*P*k18.*Q+k05.*B*k08*k12*k02*k04*k14*k16.*P*k18.*Q+k08*k12*k09.*B*k14*k01.*A*k04*k19*k21+k08*k12*k09.*B*k01.*A*k04*k19*k21*k15+k08*k12*k09.*B*k01.*A*k04*k19*k21*k17+k08*k12*k14*k16.*P*k20.*Q*k02*k04*k21+k08*k12*k14*k18.*Q*k22.*P*k02*k04*k19+k08*k12*k09.*B*k14*k01.*A*k18.*Q*k04*k19+k08*k12*k09.*B*k01.*A*k04*k19*k15*k18.*Q+k08*k12*k14*k16.*P*k18.*Q*k20.*Q*k02*k04+k08*k12*k09.*B*k14*k01.*A*k16.*P*k04*k21+k08*k12*k09.*B*k01.*A*k04*k21*k17*k16.*P+k08*k12*k14*k16.*P*k18.*Q*k22.*P*k02*k04+k08*k12*k09.*B*k14*k01.*A*k16.*P*k18.*Q*k04+k05.*B*k08*k12*k04*k19*k21*k09.*B*k14+k05.*B*k08*k12*k04*k19*k21*k09.*B*k15+k05.*B*k08*k12*k04*k19*k21*k09.*B*k17+k05.*B*k08*k04*k19*k21*k15*k13*k09.*B+k05.*B*k08*k04*k19*k21*k17*k13*k09.*B+k05.*B*k08*k12*k04*k19*k09.*B*k14*k18.*Q+k05.*B*k08*k12*k04*k19*k09.*B*k15*k18.*Q+k05.*B*k08*k04*k19*k15*k13*k18.*Q*k09.*B+k05.*B*k08*k12*k04*k21*k09.*B*k14*k16.*P+k05.*B*k08*k12*k04*k21*k09.*B*k17*k16.*P+k05.*B*k08*k04*k21*k17*k13*k16.*P*k09.*B+k05.*B*k08*k12*k04*k09.*B*k14*k16.*P*k18.*Q+k08*k12*k09.*B*k14*k16.*P*k20.*Q*k04*k21+k08*k12*k09.*B*k14*k18.*Q*k22.*P*k04*k19+k08*k12*k09.*B*k14*k16.*P*k18.*Q*k20.*Q*k04+k08*k12*k09.*B*k14*k16.*P*k18.*Q*k22.*P*k04;
E5 = k07.*I*k05.*B*k02*k04*k19*k21*k10*k14+k07.*I*k05.*B*k02*k04*k19*k21*k10*k15+k07.*I*k05.*B*k02*k04*k19*k21*k10*k17+k07.*I*k05.*B*k12*k02*k04*k19*k21*k14+k07.*I*k05.*B*k12*k02*k04*k19*k21*k15+k07.*I*k05.*B*k12*k02*k04*k19*k21*k17+k07.*I*k05.*B*k02*k04*k19*k21*k15*k13+k07.*I*k05.*B*k02*k04*k19*k21*k17*k13+k07.*I*k05.*B*k02*k04*k19*k10*k14*k18.*Q+k07.*I*k05.*B*k02*k04*k19*k10*k15*k18.*Q+k07.*I*k05.*B*k12*k02*k04*k19*k14*k18.*Q+k07.*I*k05.*B*k12*k02*k04*k19*k15*k18.*Q+k07.*I*k05.*B*k02*k04*k19*k15*k13*k18.*Q+k07.*I*k05.*B*k02*k04*k21*k10*k14*k16.*P+k07.*I*k05.*B*k02*k04*k21*k10*k17*k16.*P+k07.*I*k05.*B*k12*k02*k04*k21*k14*k16.*P+k07.*I*k05.*B*k12*k02*k04*k21*k17*k16.*P+k07.*I*k05.*B*k02*k04*k21*k17*k13*k16.*P+k07.*I*k05.*B*k02*k04*k10*k14*k16.*P*k18.*Q+k07.*I*k05.*B*k12*k02*k04*k14*k16.*P*k18.*Q+k07.*I*k12*k09.*B*k14*k01.*A*k04*k19*k21+k07.*I*k12*k09.*B*k01.*A*k04*k19*k21*k15+k07.*I*k12*k09.*B*k01.*A*k04*k19*k21*k17+k07.*I*k12*k14*k16.*P*k20.*Q*k02*k04*k21+k07.*I*k12*k14*k18.*Q*k22.*P*k02*k04*k19+k07.*I*k12*k09.*B*k14*k01.*A*k18.*Q*k04*k19+k07.*I*k12*k09.*B*k01.*A*k04*k19*k15*k18.*Q+k07.*I*k12*k14*k16.*P*k18.*Q*k20.*Q*k02*k04+k07.*I*k12*k09.*B*k14*k01.*A*k16.*P*k04*k21+k07.*I*k12*k09.*B*k01.*A*k04*k21*k17*k16.*P+k07.*I*k12*k14*k16.*P*k18.*Q*k22.*P*k02*k04+k07.*I*k12*k09.*B*k14*k01.*A*k16.*P*k18.*Q*k04+k07.*I*k05.*B*k12*k04*k19*k21*k09.*B*k14+k07.*I*k05.*B*k12*k04*k19*k21*k09.*B*k15+k07.*I*k05.*B*k12*k04*k19*k21*k09.*B*k17+k07.*I*k05.*B*k04*k19*k21*k15*k13*k09.*B+k07.*I*k05.*B*k04*k19*k21*k17*k13*k09.*B+k07.*I*k05.*B*k12*k04*k19*k09.*B*k14*k18.*Q+k07.*I*k05.*B*k12*k04*k19*k09.*B*k15*k18.*Q+k07.*I*k05.*B*k04*k19*k15*k13*k18.*Q*k09.*B+k07.*I*k05.*B*k12*k04*k21*k09.*B*k14*k16.*P+k07.*I*k05.*B*k12*k04*k21*k09.*B*k17*k16.*P+k07.*I*k05.*B*k04*k21*k17*k13*k16.*P*k09.*B+k07.*I*k05.*B*k12*k04*k09.*B*k14*k16.*P*k18.*Q+k07.*I*k12*k09.*B*k14*k16.*P*k20.*Q*k04*k21+k07.*I*k12*k09.*B*k14*k18.*Q*k22.*P*k04*k19+k07.*I*k12*k09.*B*k14*k16.*P*k18.*Q*k20.*Q*k04+k07.*I*k12*k09.*B*k14*k16.*P*k18.*Q*k22.*P*k04;
E6 = k09.*B*k14*k01.*A*k04*k06*k19*k21*k08+k09.*B*k01.*A*k04*k06*k19*k21*k08*k15+k09.*B*k01.*A*k04*k06*k19*k21*k08*k17+k11.*A*k14*k05.*B*k08*k02*k04*k19*k21+k11.*A*k05.*B*k08*k02*k04*k19*k21*k15+k11.*A*k05.*B*k08*k02*k04*k19*k21*k17+k14*k16.*P*k20.*Q*k02*k04*k06*k21*k08+k14*k18.*Q*k22.*P*k02*k04*k06*k19*k08+k09.*B*k14*k01.*A*k18.*Q*k04*k06*k19*k08+k09.*B*k01.*A*k04*k06*k19*k08*k15*k18.*Q+k11.*A*k14*k05.*B*k08*k18.*Q*k02*k04*k19+k11.*A*k05.*B*k08*k02*k04*k19*k15*k18.*Q+k14*k16.*P*k18.*Q*k20.*Q*k02*k04*k06*k08+k09.*B*k14*k01.*A*k16.*P*k04*k06*k21*k08+k09.*B*k01.*A*k04*k06*k21*k08*k17*k16.*P+k11.*A*k14*k05.*B*k08*k16.*P*k02*k04*k21+k11.*A*k05.*B*k08*k02*k04*k21*k17*k16.*P+k14*k16.*P*k18.*Q*k22.*P*k02*k04*k06*k08+k09.*B*k14*k01.*A*k16.*P*k18.*Q*k04*k06*k08+k11.*A*k14*k05.*B*k08*k16.*P*k18.*Q*k02*k04+k09.*B*k11.*A*k14*k01.*A*k08*k04*k19*k21+k09.*B*k11.*A*k01.*A*k08*k04*k19*k21*k15+k09.*B*k11.*A*k01.*A*k08*k04*k19*k21*k17+k11.*A*k14*k08*k16.*P*k20.*Q*k02*k04*k21+k11.*A*k14*k08*k18.*Q*k22.*P*k02*k04*k19+k09.*B*k11.*A*k14*k01.*A*k08*k18.*Q*k04*k19+k09.*B*k11.*A*k01.*A*k08*k04*k19*k15*k18.*Q+k11.*A*k14*k08*k16.*P*k18.*Q*k20.*Q*k02*k04+k09.*B*k11.*A*k14*k01.*A*k08*k16.*P*k04*k21+k09.*B*k11.*A*k01.*A*k08*k04*k21*k17*k16.*P+k11.*A*k14*k08*k16.*P*k18.*Q*k22.*P*k02*k04+k09.*B*k11.*A*k14*k01.*A*k08*k16.*P*k18.*Q*k04+k09.*B*k11.*A*k14*k05.*B*k08*k04*k19*k21+k09.*B*k11.*A*k05.*B*k08*k04*k19*k21*k15+k09.*B*k11.*A*k05.*B*k08*k04*k19*k21*k17+k09.*B*k14*k16.*P*k20.*Q*k04*k06*k21*k08+k09.*B*k14*k18.*Q*k22.*P*k04*k06*k19*k08+k09.*B*k11.*A*k14*k05.*B*k08*k18.*Q*k04*k19+k09.*B*k11.*A*k05.*B*k08*k04*k19*k15*k18.*Q+k09.*B*k14*k16.*P*k18.*Q*k20.*Q*k04*k06*k08+k09.*B*k11.*A*k14*k05.*B*k08*k16.*P*k04*k21+k09.*B*k11.*A*k05.*B*k08*k04*k21*k17*k16.*P+k09.*B*k14*k16.*P*k18.*Q*k22.*P*k04*k06*k08+k09.*B*k11.*A*k14*k05.*B*k08*k16.*P*k18.*Q*k04+k09.*B*k11.*A*k14*k08*k16.*P*k20.*Q*k04*k21+k09.*B*k11.*A*k14*k08*k18.*Q*k22.*P*k04*k19+k09.*B*k11.*A*k14*k08*k16.*P*k18.*Q*k20.*Q*k04+k09.*B*k11.*A*k14*k08*k16.*P*k18.*Q*k22.*P*k04;
E7 = k13*k09.*B*k01.*A*k04*k06*k19*k21*k08+k16.*P*k20.*Q*k02*k04*k06*k21*k10*k08+k18.*Q*k22.*P*k02*k04*k06*k19*k10*k08+k13*k11.*A*k05.*B*k08*k02*k04*k19*k21+k16.*P*k20.*Q*k02*k04*k06*k21*k08*k12+k18.*Q*k22.*P*k02*k04*k06*k19*k08*k12+k13*k16.*P*k20.*Q*k02*k04*k06*k21*k08+k13*k18.*Q*k22.*P*k02*k04*k06*k19*k08+k13*k18.*Q*k09.*B*k01.*A*k04*k06*k19*k08+k16.*P*k18.*Q*k20.*Q*k02*k04*k06*k10*k08+k13*k18.*Q*k11.*A*k05.*B*k08*k02*k04*k19+k16.*P*k18.*Q*k20.*Q*k02*k04*k06*k08*k12+k13*k16.*P*k18.*Q*k20.*Q*k02*k04*k06*k08+k13*k16.*P*k09.*B*k01.*A*k04*k06*k21*k08+k16.*P*k18.*Q*k22.*P*k02*k04*k06*k10*k08+k13*k16.*P*k11.*A*k05.*B*k08*k02*k04*k21+k16.*P*k18.*Q*k22.*P*k02*k04*k06*k08*k12+k13*k16.*P*k18.*Q*k22.*P*k02*k04*k06*k08+k13*k16.*P*k18.*Q*k09.*B*k01.*A*k04*k06*k08+k13*k16.*P*k18.*Q*k11.*A*k05.*B*k08*k02*k04+k13*k09.*B*k11.*A*k01.*A*k08*k04*k19*k21+k16.*P*k20.*Q*k02*k04*k21*k10*k11.*A*k08+k18.*Q*k22.*P*k02*k04*k19*k10*k11.*A*k08+k13*k16.*P*k11.*A*k20.*Q*k08*k02*k04*k21+k13*k18.*Q*k11.*A*k22.*P*k08*k02*k04*k19+k13*k18.*Q*k09.*B*k11.*A*k01.*A*k08*k04*k19+k16.*P*k18.*Q*k20.*Q*k02*k04*k10*k11.*A*k08+k13*k16.*P*k18.*Q*k11.*A*k20.*Q*k08*k02*k04+k13*k16.*P*k09.*B*k11.*A*k01.*A*k08*k04*k21+k16.*P*k18.*Q*k22.*P*k02*k04*k10*k11.*A*k08+k13*k16.*P*k18.*Q*k11.*A*k22.*P*k08*k02*k04+k13*k16.*P*k18.*Q*k09.*B*k11.*A*k01.*A*k08*k04+k13*k09.*B*k11.*A*k05.*B*k08*k04*k19*k21+k16.*P*k20.*Q*k04*k06*k21*k08*k12*k09.*B+k18.*Q*k22.*P*k04*k06*k19*k08*k12*k09.*B+k13*k16.*P*k09.*B*k20.*Q*k04*k06*k21*k08+k13*k18.*Q*k09.*B*k22.*P*k04*k06*k19*k08+k13*k18.*Q*k09.*B*k11.*A*k05.*B*k08*k04*k19+k16.*P*k18.*Q*k20.*Q*k04*k06*k08*k12*k09.*B+k13*k16.*P*k18.*Q*k09.*B*k20.*Q*k04*k06*k08+k13*k16.*P*k09.*B*k11.*A*k05.*B*k08*k04*k21+k16.*P*k18.*Q*k22.*P*k04*k06*k08*k12*k09.*B+k13*k16.*P*k18.*Q*k09.*B*k22.*P*k04*k06*k08+k13*k16.*P*k18.*Q*k09.*B*k11.*A*k05.*B*k08*k04+k13*k16.*P*k09.*B*k11.*A*k20.*Q*k08*k04*k21+k13*k18.*Q*k09.*B*k11.*A*k22.*P*k08*k04*k19+k13*k16.*P*k18.*Q*k09.*B*k11.*A*k20.*Q*k08*k04+k13*k16.*P*k18.*Q*k09.*B*k11.*A*k22.*P*k08*k04;
E8 = k20.*Q*k02*k04*k06*k21*k10*k08*k14+k20.*Q*k15*k02*k04*k06*k21*k10*k08+k20.*Q*k02*k04*k06*k21*k10*k08*k17+k20.*Q*k02*k04*k06*k21*k08*k12*k14+k20.*Q*k15*k02*k04*k06*k21*k08*k12+k20.*Q*k02*k04*k06*k21*k08*k12*k17+k20.*Q*k15*k02*k04*k06*k21*k13*k08+k20.*Q*k02*k04*k06*k21*k08*k17*k13+k20.*Q*k02*k04*k06*k10*k08*k14*k18.*Q+k20.*Q*k15*k02*k04*k06*k18.*Q*k10*k08+k20.*Q*k02*k04*k06*k08*k12*k14*k18.*Q+k20.*Q*k15*k02*k04*k06*k18.*Q*k08*k12+k20.*Q*k15*k02*k04*k06*k13*k18.*Q*k08+k15*k13*k09.*B*k01.*A*k04*k06*k21*k08+k15*k18.*Q*k22.*P*k02*k04*k06*k10*k08+k15*k13*k11.*A*k05.*B*k08*k02*k04*k21+k15*k18.*Q*k22.*P*k02*k04*k06*k08*k12+k15*k13*k18.*Q*k22.*P*k02*k04*k06*k08+k15*k13*k18.*Q*k09.*B*k01.*A*k04*k06*k08+k15*k13*k18.*Q*k11.*A*k05.*B*k08*k02*k04+k20.*Q*k02*k04*k21*k10*k11.*A*k14*k08+k20.*Q*k15*k02*k04*k21*k10*k11.*A*k08+k20.*Q*k02*k04*k21*k10*k17*k11.*A*k08+k20.*Q*k15*k02*k04*k21*k13*k11.*A*k08+k20.*Q*k02*k04*k21*k17*k13*k11.*A*k08+k20.*Q*k02*k04*k10*k11.*A*k14*k08*k18.*Q+k20.*Q*k15*k02*k04*k18.*Q*k10*k11.*A*k08+k20.*Q*k15*k02*k04*k13*k18.*Q*k11.*A*k08+k15*k13*k09.*B*k11.*A*k01.*A*k08*k04*k21+k15*k18.*Q*k22.*P*k02*k04*k10*k11.*A*k08+k15*k13*k18.*Q*k11.*A*k22.*P*k08*k02*k04+k15*k13*k18.*Q*k09.*B*k11.*A*k01.*A*k08*k04+k20.*Q*k04*k06*k21*k08*k12*k09.*B*k14+k20.*Q*k15*k04*k06*k21*k08*k12*k09.*B+k20.*Q*k04*k06*k21*k08*k12*k17*k09.*B+k20.*Q*k15*k04*k06*k21*k13*k08*k09.*B+k20.*Q*k04*k06*k21*k08*k17*k13*k09.*B+k20.*Q*k04*k06*k08*k12*k09.*B*k14*k18.*Q+k20.*Q*k15*k04*k06*k18.*Q*k08*k12*k09.*B+k20.*Q*k15*k04*k06*k13*k18.*Q*k08*k09.*B+k15*k13*k09.*B*k11.*A*k05.*B*k08*k04*k21+k15*k18.*Q*k22.*P*k04*k06*k08*k12*k09.*B+k15*k13*k18.*Q*k09.*B*k22.*P*k04*k06*k08+k15*k13*k18.*Q*k09.*B*k11.*A*k05.*B*k08*k04+k20.*Q*k15*k04*k21*k13*k09.*B*k11.*A*k08+k20.*Q*k04*k21*k17*k13*k09.*B*k11.*A*k08+k20.*Q*k15*k04*k13*k18.*Q*k09.*B*k11.*A*k08+k15*k13*k18.*Q*k09.*B*k11.*A*k22.*P*k08*k04;
E9 = k22.*P*k02*k04*k06*k19*k10*k08*k14+k22.*P*k02*k04*k06*k19*k10*k08*k15+k22.*P*k17*k02*k04*k06*k19*k10*k08+k22.*P*k02*k04*k06*k19*k08*k12*k14+k22.*P*k02*k04*k06*k19*k08*k12*k15+k22.*P*k17*k02*k04*k06*k19*k08*k12+k22.*P*k02*k04*k06*k19*k08*k15*k13+k22.*P*k17*k02*k04*k06*k19*k13*k08+k17*k13*k09.*B*k01.*A*k04*k06*k19*k08+k17*k16.*P*k20.*Q*k02*k04*k06*k10*k08+k17*k13*k11.*A*k05.*B*k08*k02*k04*k19+k17*k16.*P*k20.*Q*k02*k04*k06*k08*k12+k17*k13*k16.*P*k20.*Q*k02*k04*k06*k08+k22.*P*k02*k04*k06*k10*k08*k14*k16.*P+k22.*P*k17*k02*k04*k06*k16.*P*k10*k08+k22.*P*k02*k04*k06*k08*k12*k14*k16.*P+k22.*P*k17*k02*k04*k06*k16.*P*k08*k12+k22.*P*k17*k02*k04*k06*k13*k16.*P*k08+k17*k13*k16.*P*k09.*B*k01.*A*k04*k06*k08+k17*k13*k16.*P*k11.*A*k05.*B*k08*k02*k04+k22.*P*k02*k04*k19*k10*k11.*A*k14*k08+k22.*P*k02*k04*k19*k10*k15*k11.*A*k08+k22.*P*k17*k02*k04*k19*k10*k11.*A*k08+k22.*P*k02*k04*k19*k15*k13*k11.*A*k08+k22.*P*k17*k02*k04*k19*k13*k11.*A*k08+k17*k13*k09.*B*k11.*A*k01.*A*k08*k04*k19+k17*k16.*P*k20.*Q*k02*k04*k10*k11.*A*k08+k17*k13*k16.*P*k11.*A*k20.*Q*k08*k02*k04+k22.*P*k02*k04*k10*k11.*A*k14*k08*k16.*P+k22.*P*k17*k02*k04*k16.*P*k10*k11.*A*k08+k22.*P*k17*k02*k04*k13*k16.*P*k11.*A*k08+k17*k13*k16.*P*k09.*B*k11.*A*k01.*A*k08*k04+k22.*P*k04*k06*k19*k08*k12*k09.*B*k14+k22.*P*k04*k06*k19*k08*k12*k15*k09.*B+k22.*P*k17*k04*k06*k19*k08*k12*k09.*B+k22.*P*k04*k06*k19*k08*k15*k13*k09.*B+k22.*P*k17*k04*k06*k19*k13*k08*k09.*B+k17*k13*k09.*B*k11.*A*k05.*B*k08*k04*k19+k17*k16.*P*k20.*Q*k04*k06*k08*k12*k09.*B+k17*k13*k16.*P*k09.*B*k20.*Q*k04*k06*k08+k22.*P*k04*k06*k08*k12*k09.*B*k14*k16.*P+k22.*P*k17*k04*k06*k16.*P*k08*k12*k09.*B+k22.*P*k17*k04*k06*k13*k16.*P*k08*k09.*B+k17*k13*k16.*P*k09.*B*k11.*A*k05.*B*k08*k04+k22.*P*k04*k19*k15*k13*k09.*B*k11.*A*k08+k22.*P*k17*k04*k19*k13*k09.*B*k11.*A*k08+k17*k13*k16.*P*k09.*B*k11.*A*k20.*Q*k08*k04+k22.*P*k17*k04*k13*k16.*P*k09.*B*k11.*A*k08;
% Denominator terms
D = E1+E2+E3+E4+E5+E6+E7+E8+E9;
% Enzyme abundances terms
E1 = E1./D;
E2 = E2./D;
E3 = E3./D;
E4 = E4./D;
E5 = E5./D;
E6 = E6./D;
E7 = E7./D;
E8 = E8./D;
E9 = E9./D;
% Reaction rate 
v = +k21.*E9-k22.*P.*E1+k15.*E7-k16.*P.*E8;