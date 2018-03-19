function [ output_args ] = export_surface_handle( h_call,event_data,surf_handle,map,filename )


c_data=get(surf_handle,'CDATA');
[X]=mat2rgb(c_data,map);


imwrite(X,filename{1});

disp(strcat('surface saved to: ',filename))
end

