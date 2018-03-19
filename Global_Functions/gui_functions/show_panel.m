function [ output_args ] = show_panel( serp,alerp,panel_title,f_1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 % find all panel handles 
q=get(f_1,'Children');
class_mat={zeros(size(q))};


for j=1:length(q)
    
    %string=class(q(j));
    class_mat(j)={class(q(j))};
    
    
end

d=find(strcmp(class_mat','matlab.ui.container.Panel'));

panels=q(d);

for k=1:length(panels)
    
   panel_title_string_compare=panels(k).Title;
   
   if (~strcmp(panel_title,panel_title_string_compare))
       
       panels(k).Visible='off';
   else
       
              panels(k).Visible='on';

   end
    
end


end

