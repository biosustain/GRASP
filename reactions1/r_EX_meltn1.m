function v = r_EX_meltn1(S,P,K)
% Mass action definition 
v = K(1)*prod(S,1)-K(2)*prod(P,1);