function [ output_args ] = refresh_defaults( h_call,event_data,h )
% this f function puts a default start stop and num stwps for each input
% keeps with the idea all 



in_variables_structure=get(h.in_variables_listbox,'UserData');


listbox_value=get(h.in_variables_listbox,'Value');



set(h.in_table_start_steps_stop,'Data',cell2mat(in_variables_structure(listbox_value,3)));

end

