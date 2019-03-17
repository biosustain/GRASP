function [f,grad] = Glycolysis_Kinetics1(x,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
v = zeros(16,26);
E = zeros(16,26);
x = [x,x(:,ones(1,25)) + diag(h*1i*ones(25,1))];
else
v = zeros(16,size(x,2));
E = zeros(16,size(x,2));
end
% Defining metabolite and enzyme species
m_atp_c = x(1,:);
m_glc_D_c = x(2,:);
m_g6p_c = x(3,:);
m_f6p_c = x(4,:);
m_fdp_c = x(5,:);
m_dhap_c = x(6,:);
m_g3p_c = x(7,:);
m_13dpg_c = x(8,:);
m_3pg_c = x(9,:);
m_2pg_c = x(10,:);
m_pep_c = x(11,:);
m_glc_D_e = x(12,:);
E(1,:) = x(13,:);
E(2,:) = x(14,:);
E(3,:) = x(15,:);
E(4,:) = x(16,:);
E(5,:) = x(17,:);
E(6,:) = x(18,:);
E(7,:) = x(19,:);
E(8,:) = x(20,:);
E(9,:) = x(21,:);
E(10,:) = x(22,:);
E(11,:) = x(23,:);
E(12,:) = x(24,:);
E(13,:) = x(25,:);
E(kinInactRxns,:) = fixedExch(:,ones(1,size(x,2)));
% Reaction rates
v(1,:) = r_HEX11([m_atp_c;m_glc_D_c;m_g6p_c;ones(1,size(x,2))],[m_g6p_c;m_f6p_c],model.rxnParams(1).kineticParams,model.rxnParams(1).KnegEff,model.rxnParams(1).L,subunits(1));
v(2,:) = r_HEX21([m_atp_c;m_glc_D_c;m_g6p_c;ones(1,size(x,2))],[m_g6p_c],model.rxnParams(2).kineticParams,model.rxnParams(2).KnegEff,model.rxnParams(2).L,subunits(2));
v(3,:) = r_PGI1([m_g6p_c;m_pep_c;m_f6p_c;m_fdp_c],model.rxnParams(3).kineticParams,model.rxnParams(3).L,subunits(3));
v(4,:) = r_PFKL1([m_f6p_c;m_atp_c;ones(1,size(x,2));m_fdp_c],[m_atp_c],[m_fdp_c],model.rxnParams(4).kineticParams,model.rxnParams(4).KnegEff,model.rxnParams(4).KposEff,model.rxnParams(4).L,subunits(4));
v(5,:) = r_PFKM1([m_f6p_c;m_atp_c;ones(1,size(x,2));m_fdp_c],[m_atp_c],[m_fdp_c],model.rxnParams(5).kineticParams,model.rxnParams(5).KnegEff,model.rxnParams(5).KposEff,model.rxnParams(5).L,subunits(5));
v(6,:) = r_PFKP1([m_f6p_c;m_atp_c;ones(1,size(x,2));m_fdp_c],[m_atp_c;m_pep_c],model.rxnParams(6).kineticParams,model.rxnParams(6).KnegEff,model.rxnParams(6).L,subunits(6));
v(7,:) = r_FBA1([m_fdp_c;m_g3p_c;m_dhap_c],model.rxnParams(7).kineticParams,model.rxnParams(7).L,subunits(7));
v(8,:) = r_GAPDH1([m_g3p_c;ones(1,size(x,2));ones(1,size(x,2));m_13dpg_c;],model.rxnParams(8).kineticParams,model.rxnParams(8).L,subunits(8));
v(9,:) = r_TPI1([m_dhap_c;m_g3p_c],model.rxnParams(9).kineticParams,model.rxnParams(9).L,subunits(9));
v(10,:) = r_PGK1([m_atp_c;m_13dpg_c;ones(1,size(x,2));m_3pg_c],model.rxnParams(10).kineticParams);
v(11,:) = r_PGM1([m_3pg_c;m_2pg_c],model.rxnParams(11).kineticParams,model.rxnParams(11).L,subunits(11));
v(12,:) = r_ENO1([m_2pg_c;m_pep_c],model.rxnParams(12).kineticParams,model.rxnParams(12).L,subunits(12));
v(13,:) = r_GLUT11(m_glc_D_e,m_glc_D_c,model.rxnParams(13).kineticParams);
v(14,:) = r_EX_atp1([],[]);
v(15,:) = r_EX_glc1([],[]);
v(16,:) = r_EX_pep1([],[]);
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