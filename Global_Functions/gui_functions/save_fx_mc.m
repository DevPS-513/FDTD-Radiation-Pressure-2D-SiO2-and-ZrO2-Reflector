function [ output_args ] = save_fx_mc( h_call,event_data,h)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 %SI UNITS

meters=1;           % [m]   are 1   
cm=meters/100;      % [cm]  are [m]/100
Teslas=1;           % [T]   Base unit is tesla, need to convert any 'gauss' units
amps=1;             % [A]   is 1
celcius=1;          % [C*]  Temperature is in Celcius
newtons=1;          % [N]   Force is in Newtons
henries=1;
Watts=1;            % [W] units of power    
mu_o=4*pi*1e-7;     % permability of free space

% CONVERSIONS
mH=henries*1000;
gauss=.0001*Teslas;                         % gauss to [T]                                         
k_orsted=1000*79.57747154594*amps/meters;   % [Kor] to [A/m]
inches=.0254*meters;                        % inches to meters
pounds=4.44822162*newtons;                  % pounds to newtons

%% LOAD COIL
% coil structure is [default value,display name, variable name]
coil_structure=get(h.coil_table,'UserData');    % get structure, third column has names
coil_table_data=get(h.coil_table,'Data');       % get user-defined gui values

% sometimes coil table is returned  as a cell, if this
% happens, convert to matrix

if iscell(coil_table_data)
coil_table_data=cell2mat(coil_table_data);        
end

% now run through every coil table value ( which was made using its
% corresponding structure) and define the variable.

[Nvar, Ncol]=size(coil_structure); % run through rows
AWG=[14:1:32];

OD_mm=[1.62814,1.45034,1.29032,1.15062,1.02362,0.91186,0.8128,0.7239,...
       0.64262,0.57404,0.51054,0.45466,0.40386,0.36068,0.32004,0.28702,0.254,0.22606,0.2032];
   
OD_insulation=[ 0.06695,0.05975,0.0534,0.0478,0.04275,0.0382,0.03425,0.0306,0.02735,0.0246,0.02205,...
    0.0197,0.01765,0.0159,0.01425,0.0128,0.01145,0.0103,0.00935]*inches; % convert to meters



for j=1:Nvar
    
    name_of_variable=char(coil_structure(j,3)); % get name from 3rd column
    variable_value=coil_table_data(j,1);        % get value from GUI table   


 eval(strcat(name_of_variable,'=',num2str(variable_value))); % var=value
 
[val,loc]=find_var(coil_structure,coil_structure(:,1),name_of_variable,3);
r_wire=interp1(AWG,OD_mm/2/1000,gauge);
OD_in=interp1(AWG,OD_insulation,gauge);
wire_buff=OD_in/2-r_wire;
eval(strcat('coil_table_values(loc)=',char(name_of_variable),';')); 



end
set(h.coil_table,'Data',coil_table_values');       % get user-defined gui values

magnet_value=get(h.magnet_type_popup,'Value');                  % get magnet name location in listbox
magnet_list_structure=get(h.magnet_type_popup,'UserData');      % get the variable name
magnet_list=magnet_list_structure(:,1);                         % get all string names
magnet_name=char(magnet_list(magnet_value));                    % get the current magnet name
num_turns_c=round(Lc./(2*(r_wire+wire_buff)));          % number of turns in coil


file_name_to_save=strcat('./saved_mat_files/',magnet_name,'_ID_',num2str(2*1000*rc,2),'_mm',...
    '_Nt_',num2str(num_turns_c),'_Ns_',num2str(N_stack),'_Nm_',num2str(N_mags),'_AWG_',num2str(gauge));

fx_mc=h.fx_mc.YData;
x_mc=h.fx_mc.XData;

% get last data in single coil-magnet output table
% assume it is recent and save stuff

currentrow=get(h.output_table,'Userdata');
gui_output_table_data=get(h.output_table,'data');
r_mat=ones(1,N_stack).*(rc+r_wire+wire_buff)+[0:1:N_stack-1]*2*(r_wire+wire_buff);

R_th=sum(resistivity*2*pi*r_mat.*num_turns_c/(pi*r_wire^2));% resistnace all stacks

Pdc=V_src^2/R_th;

save(file_name_to_save,'fx_mc','x_mc','Pdc','V_src','R_th','N_mags');
disp(strcat('saved to: ',file_name_to_save))
end

