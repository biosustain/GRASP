function v = testMassAction1(SC,S,PC,P,K)
% Mass action definition 
v = K(1)*prod(S.^SC, 1)-K(2)*prod(P.^PC, 1);
