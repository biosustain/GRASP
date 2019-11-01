function [f,grad] = toy_model_Kinetics1(x,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
v = zeros(14,8);
E = zeros(14,8);
x = [x,x(:,ones(1,7)) + diag(h*1i*ones(7,1))];
else
v = zeros(14,size(x,2));
E = zeros(14,size(x,2));
end
% Defining metabolite and enzyme species
m_m5 = x(1,:);
m_m6 = x(2,:);
m_m7 = x(3,:);
m_m8 = x(4,:);
m_m9 = x(5,:);
m_m10 = x(6,:);
m_m11 = x(7,:);
E(kinInactRxns,:) = fixedExch(:,ones(1,size(x,2)));
% Reaction rates
v(1,:) = r_r11([ones(1,size(x,2));m_m6;m_m6;m_m5;ones(1,size(x,2))],model.rxnParams(1).kineticParams);
v(2,:) = r_r21([m_m5;m_m6;m_m7;m_m10],model.rxnParams(2).kineticParams);
v(3,:) = r_r31([ones(1,size(x,2));m_m7;ones(1,size(x,2));m_m10;m_m9;m_m8;m_m11;ones(1,size(x,2));ones(1,size(x,2))],model.rxnParams(3).kineticParams);
v(4,:) = r_r4_11([ones(1,size(x,2));m_m8;m_m9;ones(1,size(x,2));],model.rxnParams(4).kineticParams);
v(5,:) = r_r4_21([ones(1,size(x,2));m_m8;m_m9;ones(1,size(x,2));],model.rxnParams(5).kineticParams);
v(6,:) = r_r51([m_m5;m_m6;m_m7;m_m10],model.rxnParams(6).kineticParams);
v(7,:) = r_r61([ones(1,size(x,2));m_m7;ones(1,size(x,2));m_m10;m_m9;m_m8;m_m11;ones(1,size(x,2));ones(1,size(x,2))],model.rxnParams(7).kineticParams);
v(8,:) = r_r71(ones(1,size(x,2)),m_m6,model.rxnParams(8).kineticParams);
v(9,:) = r_r81(m_m6,ones(1,size(x,2)),model.rxnParams(9).kineticParams);
v(10,:) = r_r91(m_m7,ones(1,size(x,2)),model.rxnParams(10).kineticParams);
v(11,:) = r_r101(m_m5,ones(1,size(x,2)),model.rxnParams(11).kineticParams);
v(12,:) = r_r111([m_m8],model.rxnParams(12).kineticParams);
v(13,:) = r_r121(m_m9,ones(1,size(x,2)),model.rxnParams(13).kineticParams);
v(14,:) = r_r131([],[]);
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