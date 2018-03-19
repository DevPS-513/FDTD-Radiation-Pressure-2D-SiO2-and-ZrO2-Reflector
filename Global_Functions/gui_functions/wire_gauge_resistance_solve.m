function [] = wire_gauge_resistance_solve(h_call,event_data,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here





inches=.0254;
AWG=[14:1:32];

OD_mm=[1.62814,1.45034,1.29032,1.15062,1.02362,0.91186,0.8128,0.7239,...
       0.64262,0.57404,0.51054,0.45466,0.40386,0.36068,0.32004,0.28702,0.254,0.22606,0.2032];
   
OD_insulation=[ 0.06695,0.05975,0.0534,0.0478,0.04275,0.0382,0.03425,0.0306,0.02735,0.0246,0.02205,...
    0.0197,0.01765,0.0159,0.01425,0.0128,0.01145,0.0103,0.00935]*inches; % convert to meters


wire_buff=(OD_insulation-OD_mm/1000)/2;


gauge_mat=zeros(size(AWG));
R_mat=zeros(size(AWG));
Ns_mat=zeros(size(AWG));
Nt_mat=zeros(size(AWG));
delta_R_mat=zeros(size(AWG));

R_sim_method_mat=zeros(size(AWG));

% get coil parameters


%% LOAD COIL
% coil structure is  [default value,display name, variable name]
coil_structure=get(h.coil_table,'UserData');    % get structure, third column has names
coil_table_data=get(h.coil_table,'Data');       % get user-defined gui values
    
% using the variable name, find its location in the structure then--
[Lc,loc]=find_var(coil_structure,coil_table_data(:,1),'Lc',3);

[rc,loc]=find_var(coil_structure,coil_table_data(:,1),'rc',3);
[rho,loc]=find_var(coil_structure,coil_table_data(:,1),'resistivity',3);
% 
% Lc=cell2mat(Lc);
% rc=cell2mat(rc);
% rho=cell2mat(rho);

N_stack_mat=[1:100];


% R for Lc=.0236


R_star=get(h.resistance_tool_editbox,'String');
R_star=str2num(R_star);

R_mat=zeros(1,length(AWG));
extra_N_mat=zeros(1,length(AWG));
delta_R_mat=zeros(1,length(AWG));
R_mat_corr=zeros(1,length(AWG));
N_stack_coil_mat=zeros(1,length(AWG));
wc_mat=zeros(1,length(AWG));

N_remove=1;








for j=1:length(AWG)

gauge=AWG(j);



OD_in=interp1(AWG,OD_insulation,gauge);
rw=interp1(AWG,OD_mm,gauge)/2/1000;
wire_buff=OD_in/2-rw;



Nt=round(Lc./OD_in);
Nt_mat(j)=Nt;
for q=1:length(N_stack_mat)
    N_stack=N_stack_mat(q);
    
%     if N_stack==45
%         
%        system('PAUSE');
%        
%     end
    
  r_mat=ones(1,N_stack).*(rc+rw+wire_buff)+[0:1:N_stack-1]*2*(rw+wire_buff);


R(q)=rho*Nt.*sum(2*pi*r_mat)/(pi*rw^2);

end
[val loc]=min(abs(R-R_star));
R_mat(j)=R(loc);
N_stack_coil_mat(j)=N_stack_mat(loc);

% if there is a big difference,
loop_R=rho*2*pi*r_mat(end)/(pi*rw^2);
delta_R=R_star-R_mat(j);
delta_R_mat(j)=delta_R;

if abs(R_star-R_mat(j))>loop_R;
extra_N_mat(j)=round(delta_R*pi*rw^2/(rho*2*pi*r_mat(end)));    
end

wc_mat(j)=N_stack_coil_mat(j)*2*(rw+wire_buff);  % height of coil

R_mat_corr(j)=R_mat(j)+rho*2*pi*extra_N_mat(j)*r_mat(end)/(pi*rw^2);


if extra_N_mat(j)<0
N_remove=extra_N_mat(j);    
end

end 

data_output_table= zeros(length(AWG),6);

data_output_table(:,1)=AWG;
data_output_table(:,2)=N_stack_coil_mat;
data_output_table(:,3)=Nt_mat;
data_output_table(:,4)=extra_N_mat;
data_output_table(:,5)=R_mat_corr;
data_output_table(:,6)=R_mat;
data_output_table(:,7)=wc_mat;


set(h.f_800,'Visible','On')
set(h.r_tool_table,'Data',data_output_table)
 




                        


end

