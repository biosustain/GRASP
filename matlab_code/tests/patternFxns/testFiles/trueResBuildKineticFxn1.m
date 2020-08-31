function [f,grad] = toy_model1_Kinetics1(x,xconst,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
xconst = xconst(:);
v = zeros(13,21);
E = zeros(13,21);
x = [x,x(:,ones(1,20)) + diag(h*1i*ones(20,1))];
xconst = [xconst,xconst(:,ones(1,20))];
else
v = zeros(13,size(x,2));
E = zeros(13,size(x,2));
end
% Defining metabolite and enzyme species
m_m5 = x(1,:);
m_m6 = x(2,:);
m_m7 = x(3,:);
m_m8 = x(4,:);
m_m9 = x(5,:);
m_m10 = x(6,:);
m_m11 = x(7,:);
E(1,:) = x(8,:);
E(2,:) = x(9,:);
E(3,:) = x(10,:);
E(4,:) = x(11,:);
E(5,:) = x(12,:);
E(6,:) = x(13,:);
E(7,:) = x(14,:);
E(8,:) = x(15,:);
E(9,:) = x(16,:);
E(10,:) = x(17,:);
E(11,:) = x(18,:);
E(12,:) = x(19,:);
E(13,:) = x(20,:);
m_m1 = xconst(1,:);
m_m2 = xconst(2,:);
m_m3 = xconst(3,:);
m_m4 = xconst(4,:);
m_m12 = xconst(5,:);
m_m13 = xconst(6,:);
m_m14 = xconst(7,:);
m_m15 = xconst(8,:);
m_m16 = xconst(9,:);
m_m17 = xconst(10,:);
m_m18 = xconst(11,:);
m_m19 = xconst(12,:);
m_m20 = xconst(13,:);
% Reaction rates
v(1,:) = r_r11([m_m3;m_m6;m_m6;m_m5;m_m14],model.rxnParams(1).kineticParams);
v(2,:) = r_r21([m_m5;m_m6;m_m7;m_m10],model.rxnParams(2).kineticParams);
v(3,:) = r_r31([m_m1;m_m7;m_m1;m_m10;m_m9;m_m8;m_m11;m_m12;m_m12],model.rxnParams(3).kineticParams);
v(4,:) = r_r41([m_m2;m_m8;m_m9;m_m13;],model.rxnParams(4).kineticParams);
v(5,:) = r_r51([m_m5;m_m6;m_m7;m_m10],model.rxnParams(5).kineticParams);
v(6,:) = r_r61([m_m1;m_m7;m_m1;m_m10;m_m9;m_m8;m_m11;m_m12;m_m12],model.rxnParams(6).kineticParams);
v(7,:) = r_r71([1*ones(1,size(x,2))],[m_m4],[1*ones(1,size(x,2))],[m_m6],model.rxnParams(7).kineticParams);
v(8,:) = r_r81([1*ones(1,size(x,2))],[m_m6],[1*ones(1,size(x,2))],[m_m16],model.rxnParams(8).kineticParams);
v(9,:) = r_r91([1*ones(1,size(x,2))],[m_m7],[1*ones(1,size(x,2))],[m_m20],model.rxnParams(9).kineticParams);
v(10,:) = r_r101([1*ones(1,size(x,2))],[m_m5],[1*ones(1,size(x,2))],[m_m15],model.rxnParams(10).kineticParams);
v(11,:) = r_r111([m_m8],model.rxnParams(11).kineticParams);
v(12,:) = r_r121([1*ones(1,size(x,2))],[m_m9],[1*ones(1,size(x,2))],[m_m19],model.rxnParams(12).kineticParams);
v(13,:) = r_r131(fixedExch(1), size(x,2));
if flag==1
% Final rates
y = sum((Sred*(E.*v)).^2);
f = real(y(1));
if (nargout>1) % gradient is required
grad = imag(y(2:end))/h;
end
else
f = E.*v;
grad = [];
end