function y = toy_model_Kinetics1_ode(y,Eref,metsRefConc,model,fixedExch,Sred,kinInactRxns,subunits,flag)
v = zeros(14,size(Eref,2));
E = zeros(14,size(Eref,2));
m_m5 = y(1,:);
m_m6 = y(2,:);
m_m7 = y(3,:);
m_m8 = y(4,:);
m_m9 = y(5,:);
m_m10 = y(6,:);
m_m11 = y(7,:);
E(1,:) = Eref(1,:);
E(2,:) = Eref(2,:);
E(3,:) = Eref(3,:);
E(4,:) = Eref(4,:);
E(5,:) = Eref(5,:);
E(6,:) = Eref(6,:);
E(7,:) = Eref(7,:);
E(8,:) = Eref(8,:);
E(9,:) = Eref(9,:);
E(10,:) = Eref(10,:);
E(11,:) = Eref(11,:);
E(12,:) = Eref(12,:);
E(13,:) = Eref(13,:);
E(14,:) = Eref(14,:);
v(1,:) = r_r11([ones(1,size(y,2));m_m6;m_m6;m_m5;ones(1,size(y,2))],model.rxnParams(1).kineticParams);
v(2,:) = r_r21([m_m5;m_m6;m_m7;m_m10],model.rxnParams(2).kineticParams);
v(3,:) = r_r31([ones(1,size(y,2));m_m7;ones(1,size(y,2));m_m10;m_m9;m_m8;m_m11;ones(1,size(y,2));ones(1,size(y,2))],model.rxnParams(3).kineticParams);
v(4,:) = r_r4_11([ones(1,size(y,2));m_m8;m_m9;ones(1,size(y,2));],model.rxnParams(4).kineticParams);
v(5,:) = r_r4_21([ones(1,size(y,2));m_m8;m_m9;ones(1,size(y,2));],model.rxnParams(5).kineticParams);
v(6,:) = r_r51([m_m5;m_m6;m_m7;m_m10],model.rxnParams(6).kineticParams);
v(7,:) = r_r61([ones(1,size(y,2));m_m7;ones(1,size(y,2));m_m10;m_m9;m_m8;m_m11;ones(1,size(y,2));ones(1,size(y,2))],model.rxnParams(7).kineticParams);
v(8,:) = r_r71([1*ones(1,size(y,2))],[ones(1,size(y,2))],[1*ones(1,size(y,2))],[m_m6],model.rxnParams(8).kineticParams);
v(9,:) = r_r81([1*ones(1,size(y,2))],[m_m6],[1*ones(1,size(y,2))],[ones(1,size(y,2))],model.rxnParams(9).kineticParams);
v(10,:) = r_r91([1*ones(1,size(y,2))],[m_m7],[1*ones(1,size(y,2))],[ones(1,size(y,2))],model.rxnParams(10).kineticParams);
v(11,:) = r_r101([1*ones(1,size(y,2))],[m_m5],[1*ones(1,size(y,2))],[ones(1,size(y,2))],model.rxnParams(11).kineticParams);
v(12,:) = r_r111([m_m8],model.rxnParams(12).kineticParams);
v(13,:) = r_r121([1*ones(1,size(y,2))],[m_m9],[1*ones(1,size(y,2))],[ones(1,size(y,2))],model.rxnParams(13).kineticParams);
v(14,:) = r_r131(fixedExch(1), size(y,2));

y = (1./(metsRefConc.*10^3)) .* (Sred*(E.*v));