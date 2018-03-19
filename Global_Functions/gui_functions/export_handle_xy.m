function [ output_args ] = export_handle_xy( h_call,event_data,h,handle,prefix,norm_on )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% LOAD COIL
% coil structure is  [default value,display name, variable name]
coil_structure=get(h.coil_table,'UserData');    % get structure, third column has names
coil_table_data=get(h.coil_table,'Data');       % get user-defined gui values


% using the variable name, find its location in the structure then--
[Lc,loc]=find_var(coil_structure,coil_table_data,'Lc',3);
[r_wire_ins,loc]=find_var(coil_structure,coil_table_data,'r_wire_ins',3);

Nt=round(Lc./(2*(r_wire_ins)));          % number of turns in coil


[Ns,loc]=find_var(coil_structure,coil_table_data,'N_stack',3);
[N_mags,loc]=find_var(coil_structure,coil_table_data,'N_mags',3);
[gauge,loc]=find_var(coil_structure,coil_table_data,'gauge',3);
[rc,loc]=find_var(coil_structure,coil_table_data,'rc',3);

%% LOAD PACK DATA
pack_structure=get(h.pack_table,'UserData'); % Get structure that has both the display name and the variable name
pack_table=get(h.pack_table,'Data');         % Get structure that has both the display name and the variable name

[N_coils,loc]=find_var(pack_structure,pack_table,'N_coils',3);    
[N_sh_mags,loc]=find_var(pack_structure,pack_table,'N_sh_mags',3);

x_data=get(handle,'XDATA');
y_data=get(handle,'YDATA');

magnet_value=get(h.magnet_type_popup,'Value');                  % get magnet name location in listbox
magnet_list_structure=get(h.magnet_type_popup,'UserData');      % get the variable name
magnet_list=magnet_list_structure(:,1);                         % get all string names
magnet_name=char(magnet_list(magnet_value));                    % get the current magnet name

% check if saturated or not

sat_on=get(h.f_21_saturation_current_enable,'Value');



filename_to_write=strcat('./exported_data/',prefix,'_',magnet_name,'_ID_',num2str(2*1000*rc,2),'_mm','_Nt_',num2str(Nt),'_Ns_',num2str(Ns),'_N_mags_',...
    num2str(N_mags),'_Nc_',num2str(N_coils),'_Nshm_',num2str(N_sh_mags),'_AWG_',num2str(gauge),'_sat_',num2str(sat_on),'.dat');

write_dat(filename_to_write,x_data,y_data);

disp(filename_to_write);

end

