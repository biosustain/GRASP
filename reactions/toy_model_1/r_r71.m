function v = r_r71(S,P,K)
% Mass action definition 
v = K(1)*prod(S,1)-K(2)*prod(P,1);