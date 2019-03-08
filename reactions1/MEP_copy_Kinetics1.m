function [f,grad] = MEP_copy_Kinetics1(x,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
v = zeros(14,36);
E = zeros(14,36);
x = [x,x(:,ones(1,35)) + diag(h*1i*ones(35,1))];
else
v = zeros(14,size(x,2));
E = zeros(14,size(x,2));
end
% Defining metabolite and enzyme species
m_g3p = x(1,:);
m_pyr = x(2,:);
m_dxp = x(3,:);
m_nadph = x(4,:);
m_nadp = x(5,:);
m_fld_ox = x(6,:);
m_fld_red = x(7,:);
m_mep = x(8,:);
m_ctp = x(9,:);
m_cmp = x(10,:);
m_cdmep = x(11,:);
m_ppi = x(12,:);
m_cdpmep = x(13,:);
m_atp = x(14,:);
m_adp = x(15,:);
m_mecpp = x(16,:);
m_hmbpp = x(17,:);
m_dmpp = x(18,:);
m_ipp = x(19,:);
m_iso = x(20,:);
m_mecpp_ex = x(21,:);
m_dxp_ex = x(22,:);
E(1,:) = x(23,:);
E(2,:) = x(24,:);
E(3,:) = x(25,:);
E(4,:) = x(26,:);
E(5,:) = x(27,:);
E(6,:) = x(28,:);
E(7,:) = x(29,:);
E(8,:) = x(30,:);
E(9,:) = x(31,:);
E(10,:) = x(32,:);
E(11,:) = x(33,:);
E(12,:) = x(34,:);
E(13,:) = x(35,:);
E(kinInactRxns,:) = fixedExch(:,ones(1,size(x,2)));
% Reaction rates
v(1,:) = r_dxs1([m_g3p;m_pyr;ones(1,size(x,2));m_dxp],model.rxnParams(1).kineticParams);
v(2,:) = r_EX_dxp1(m_dxp,m_dxp_ex,model.rxnParams(2).kineticParams);
v(3,:) = r_dxr1([m_dxp;m_nadph;m_dxp;m_nadp;m_mep],model.rxnParams(3).kineticParams);
v(4,:) = r_ispD1([m_mep;m_ctp;m_cdmep;m_ppi],model.rxnParams(4).kineticParams);
v(5,:) = r_ispE1([m_cdmep;m_atp;m_cdpmep;m_adp],model.rxnParams(5).kineticParams);
v(6,:) = r_ispF1([m_cdpmep;m_cmp;m_mecpp],[m_mep],model.rxnParams(6).kineticParams,model.rxnParams(6).KposEff,model.rxnParams(6).L,subunits(6));
v(7,:) = r_EX_mecpp1([m_mecpp;m_mecpp_ex],model.rxnParams(7).kineticParams);
v(8,:) = r_ispG1([m_fld_red;m_mecpp;m_fld_ox;m_hmbpp],model.rxnParams(8).kineticParams);
v(9,:) = r_fpr1(m_nadph,m_nadp,model.rxnParams(9).kineticParams);
v(10,:) = r_ispH1([m_fld_red;m_hmbpp;m_fld_ox;m_dmpp;m_ipp],model.rxnParams(10).kineticParams);
v(11,:) = r_idi1([m_ipp;m_dmpp],model.rxnParams(11).kineticParams);
v(12,:) = r_ispS1([m_dmpp;m_ppi;m_iso],model.rxnParams(12).kineticParams);
v(13,:) = r_EX_isoprene1([m_iso],model.rxnParams(13).kineticParams);
v(14,:) = r_biomass_drain1([],[]);
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