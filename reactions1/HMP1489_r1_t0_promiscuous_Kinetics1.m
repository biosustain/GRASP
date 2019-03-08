function [f,grad] = HMP1489_r1_t0_promiscuous_Kinetics1(x,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
v = zeros(10,18);
E = zeros(10,18);
x = [x,x(:,ones(1,17)) + diag(h*1i*ones(17,1))];
else
v = zeros(10,size(x,2));
E = zeros(10,size(x,2));
end
% Defining metabolite and enzyme species
m_fivehtp_c = x(1,:);
m_trp_c = x(2,:);
m_srtn_c = x(3,:);
m_nactsertn_c = x(4,:);
m_meltn_c = x(5,:);
m_tryptm_c = x(6,:);
m_nactryptm_c = x(7,:);
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
% Reaction rates
v(1,:) = r_TPH1([ones(1,size(x,2));m_trp_c;m_trp_c;m_fivehtp_c;ones(1,size(x,2))],model.rxnParams(1).kineticParams);
v(2,:) = r_DDC1([m_fivehtp_c;m_trp_c;m_srtn_c;m_tryptm_c],model.rxnParams(2).kineticParams);
v(3,:) = r_AANAT1([ones(1,size(x,2));m_srtn_c;ones(1,size(x,2));m_tryptm_c;m_meltn_c;m_nactsertn_c;m_nactryptm_c;ones(1,size(x,2));ones(1,size(x,2))],model.rxnParams(3).kineticParams);
v(4,:) = r_ASMT1([ones(1,size(x,2));m_nactsertn_c;m_srtn_c;m_meltn_c;ones(1,size(x,2));],model.rxnParams(4).kineticParams);
v(5,:) = r_DDC_tryptm1([m_fivehtp_c;m_trp_c;m_srtn_c;m_tryptm_c],model.rxnParams(5).kineticParams);
v(6,:) = r_AANAT_tryptm1([ones(1,size(x,2));m_srtn_c;ones(1,size(x,2));m_tryptm_c;m_meltn_c;m_nactsertn_c;m_nactryptm_c;ones(1,size(x,2));ones(1,size(x,2))],model.rxnParams(6).kineticParams);
v(7,:) = r_IN_trp1(ones(1,size(x,2)),m_trp_c,model.rxnParams(7).kineticParams);
v(8,:) = r_EX_trp1(m_trp_c,ones(1,size(x,2)),model.rxnParams(8).kineticParams);
v(9,:) = r_EX_meltn1(m_meltn_c,ones(1,size(x,2)),model.rxnParams(9).kineticParams);
v(10,:) = r_EX_nactryptm1(m_nactryptm_c,ones(1,size(x,2)),model.rxnParams(10).kineticParams);
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