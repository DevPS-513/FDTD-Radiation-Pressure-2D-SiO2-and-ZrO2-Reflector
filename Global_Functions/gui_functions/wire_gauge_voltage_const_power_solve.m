function [] = wire_gauge_voltage_const_power_solve(h_call,event_data,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here





inches=.0254;
AWG=[14:1:32];

OD_mm=[1.62814,1.45034,1.29032,1.15062,1.02362,0.91186,0.8128,0.7239,...
       0.64262,0.57404,0.51054,0.45466,0.40386,0.36068,0.32004,0.28702,0.254,0.22606,0.2032];
   
OD_insulation=[ 0.06695,0.05975,0.0534,0.0478,0.04275,0.0382,0.03425,0.0306,0.02735,0.0246,0.02205,...
    0.0197,0.01765,0.0159,0.01425,0.0128,0.01145,0.0103,0.00935]*inches; % convert to meters





%% LOAD COIL
% coil structure is  [default value,display name, variable name]
coil_structure=get(h.coil_table,'UserData');    % get structure, third column has names
coil_table_data=get(h.coil_table,'Data');       % get user-defined gui values
    
% using the variable name, find its location in the structure then--
[Lc,loc]=find_var(coil_structure,coil_table_data(:,1),'Lc',3);
[rw_norm,loc]=find_var(coil_structure,coil_table_data(:,1),'r_wire_ins',3);
[N_stack_norm,loc]=find_var(coil_structure,coil_table_data(:,1),'N_stack',3);
[P_star,loc]=find_var(coil_structure,coil_table_data(:,1),'P_star',3);

[rc,loc]=find_var(coil_structure,coil_table_data(:,1),'rc',3);
[rho,loc]=find_var(coil_structure,coil_table_data(:,1),'resistivity',3);
% 
% Lc=cell2mat(Lc);
% rc=cell2mat(rc);
% rho=cell2mat(rho);


wc=N_stack_norm*2*interp1(AWG,OD_insulation,22)/2;

% R for Lc=.0236





R_mat=zeros(1,length(AWG));
V_mat=zeros(1,length(AWG));
I_mat=zeros(1,length(AWG));
N_stack_mat=zeros(1,length(AWG));

wc_mat=zeros(size(AWG));
Lc_mat=zeros(size(AWG));

Nt_mat=zeros(size(AWG));





for j=1:length(AWG)

gauge=AWG(j);

OD_in=interp1(AWG,OD_insulation,gauge);
rw=interp1(AWG,OD_mm,gauge)/2/1000;
wire_buff=OD_in/2-rw;

Nt=round(Lc./(2*(rw+wire_buff)));
N_stack=round(wc./(2*(rw+wire_buff)));

r_mat=ones(1,N_stack).*(rc+rw+wire_buff)+[0:1:N_stack-1]*2*(rw+wire_buff);

R_th=rho*Nt*sum(2*pi*r_mat)/(pi*rw^2);

R_mat(j)=R_th;
I_mat(j)=sqrt(P_star/R_th);
V_mat(j)=I_mat(j)*R_mat(j);

Nt_mat(j)=Nt;
N_stack_mat(j)=N_stack;
wc_mat(j)=N_stack*2*(rw+wire_buff);
Lc_mat(j)=Nt*2*(rw+wire_buff);

end 



f_800=figure(800);
data_output_table= zeros(length(AWG),8);

data_output_table(:,1)=AWG;
data_output_table(:,2)=Nt_mat;
data_output_table(:,3)=N_stack_mat;
data_output_table(:,4)=R_mat;
data_output_table(:,5)=I_mat;
data_output_table(:,6)=V_mat;
data_output_table(:,7)=wc_mat;
data_output_table(:,8)=Lc_mat;



h.coil_table = uitable('parent',f_800,'Data',data_output_table,...
    'ColumnName',{'AWG','Nt','Ns','R','I','V','wc','Lc'},...
    'units','normalized',...
    'position',[0,0,1,1],...
    'ColumnEditable',[true true]);
                        


end

