function allostericTest
Vmax = 10;
Keq  = 2;
Ka   = 2.5;
Kb   = 2;

% Try with simple MM kinetics
A = linspace(1,1e2,50);
B = 2;
v = (Vmax/Ka)*(A - B/Keq)./(1 + A/Ka + B/Kb);
plot(A,v)
hold on

% Make v allosteric
n = 2;
Lo = 1e-3;
v = v.*(Lo*(1 + A/Ka).^-n + 1).^-1;
plot(A,v,'-r')
hold off