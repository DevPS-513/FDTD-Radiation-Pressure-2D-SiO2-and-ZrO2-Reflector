
function [N2X,cx1,cx2,cy1,cy2,box_1x_mat,box_1y_mat,box_2x_mat,box_2y_mat]=create_bragg_grating_v2(L_x,L_y,c_x,c_y,er_1,er_2,L_1,L_2,N2X,dx,dy)
% This function creates a bragg grating of total width and height L_x by
% L_y and centered at c_x and c_y on the N2X grid
% has inputs er_1,L_1,er_2,L_2 being the first and second layer
% premitiivity and length, outputs back the geomretry N2X, the four corners
% of the grating as cx1,cx2,cy1, and cy2, and box_1x_mat and box_1y_mat 
% as corner boxes defining layer 1, i.e 5 points per box
% (x1,y1)--(x1,y2)--(x2,y2)--(x2,y1)--(x1,y1)


r_nx=round(L_x/dx/2);
r_ny=round(L_y/dy/2);

m_y=round(c_y/dy);
m_x=round(c_x/dx);

% n as a function of r
cx1=m_x-r_nx;
cx2=m_x+r_nx;

cy1=m_y-r_ny;
cy2=m_y+r_ny;

n_1=round(L_2/dx);
n_2=round(L_1/dx);


% index of metal+dielectric
n_cell=round((L_1+L_2)/dx)+1;

make_num=round(L_x/(L_1+L_2));

% boces to store the four corners of each layer
box_1x_mat=zeros(1,make_num*5);
box_1y_mat=zeros(1,make_num*5);
box_2x_mat=zeros(1,make_num*5);
box_2y_mat=zeros(1,make_num*5);
box_counter=1;

 N2X(cx1:cx2,cy1:cy2)=er_2;
       
        for make_j=1:make_num
            
            
for i_x=cy1:cy2    
       
        
            n_start=cx1+(make_j-1)*n_cell;
            no=n_start;
            n1=n_start+n_2;
            n2=n1+1;
            n3=n2+n_1;     
            
            N2X(no:n1,i_x)=er_1;
       
%             figure
%             plot([cy1:cy2]*dx*1e9,N2X(i_x,cy1:cy2))
%             hold on
%             line(n_start*dx*1E9*[1 1],[0 er_core],'color','r')
%             line((n_start+n_m)*dx*1E9*[1 1],[0 er_core],'color','g')
%             line((n_m+1+n_d)*dx*1E9*[1 1],[0 er_core],'color','black')     
               end   
        
         box_1x_mat(1,box_counter:box_counter+4)=[no no n1 n1 no];
         box_1y_mat(1,box_counter:box_counter+4)=[cy1 cy2 cy2 cy1 cy1];
         
         box_2x_mat(1,box_counter:box_counter+4)=[n2 n2 n3 n3 n2];
         box_2y_mat(1,box_counter:box_counter+4)=[cy1 cy2 cy2 cy1 cy1];
box_counter=box_counter+5;
     
      
     
end
  cx2=n3; % overwrite cx2 with last x_index to be written 