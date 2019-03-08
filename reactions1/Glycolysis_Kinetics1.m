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
