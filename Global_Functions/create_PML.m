function [ er,mr,sigma_e,sigma_m ] = create_PML(sigma_1,sigma_PML_max,er,mr,NPML_x,NPML_y,n_os_x,n_os_y )
% This functions defines a sigma PML
c=299792458;
mu_o=4*pi*10^-7;
eps_o=(1/(c*c*mu_o));

[Nx,Ny]=size(er);
sigma_e=zeros(Nx,Ny);
sigma_m=zeros(Nx,Ny);

% MATCH EPS and MU
eps_1_left=er(NPML_x+n_os_x,:);        	% Sample Surrounding Medium, call this medium 1     
mu_1_left=mr(NPML_x+n_os_x,:);

eps_1_right=er(Nx-NPML_x+1-n_os_x,:);                                     % Sample Surrounding Medium, call this medium 1     
mu_1_right=mr(Nx-NPML_x+1-n_os_x,:);

eps_1_bottom=er(:,NPML_y+n_os_y);                                         % Sample Surrounding Medium, call this medium 1     
mu_1_bottom=mr(:,NPML_y+n_os_y,:);

eps_1_top=er(:,Ny-NPML_y+1-n_os_y);                                       % Sample Surrounding Medium, call this medium 1     
mu_1_top=mr(:,Ny-NPML_y+1-n_os_y);

eps_2_left=(eps_1_left);                                                % Solve for Medium 2 in PML, assuming eps2=eps1
mu_2_left=(eps_2_left./eps_1_left).*mu_1_left;                          % This will only work if mr=1 for surrounding materials

eps_2_right=(eps_1_right);                                              % Solve for Medium 2 in PML, assuming eps2=eps1
mu_2_right=(eps_2_right./eps_1_right).*mu_1_right;                      % This will only work if mr=1 for surrounding materials

eps_2_top=(eps_1_top);                                                  % Solve for Medium 2 in PML, assuming eps2=eps1
mu_2_top=(eps_2_top./eps_1_top).*mu_1_top;                              % This will only work if mr=1 for surrounding materials

eps_2_bottom=(eps_1_bottom);                                            % Solve for Medium 2 in PML, assuming eps2=eps1
mu_2_bottom=(eps_2_bottom./eps_1_bottom).*mu_1_bottom;                  % This will only work if mr=1 for surrounding materials

    
    
m_fac_left=(mu_2_left./(eps_2_left)).*(mu_o/eps_o);             % sigma_m=m_fac*sigma_e
m_fac_right=(mu_2_right./(eps_2_right)).*(mu_o/eps_o);          % sigma_m=m_fac*sigma_e

m_fac_top=(mu_2_top./(eps_2_top)).*(mu_o/eps_o);                % sigma_m=m_fac*sigma_e
m_fac_bottom=(mu_2_bottom./(eps_2_bottom)).*(mu_o/eps_o);       % sigma_m=m_fac*sigma_e

    
    for j_d=1:Ny
  % Make left PML region equal to permitivity just to the right
    er(1:NPML_x+1+n_os_x,j_d)=er(NPML_x+1+n_os_x,j_d);
  % Make right PML region equal to permmitivity right to the left 
    er((Nx-NPML_x-1-(n_os_x)):Nx,j_d)=er(Nx-NPML_x-1-n_os_x,j_d);
  % Make left PML mr region equal to mr just to the right
    mr(1:NPML_x+1+n_os_x,j_d)=mr(NPML_x+1+n_os_x,j_d);
  % Make right PML mr region equal to mr just to the left
    mr((Nx-NPML_x-1-(n_os_x)):Nx,j_d)=mr(Nx-NPML_x-1-n_os_x,j_d);        
    end
    
  
    
    for i_d=1:Nx
% bottom er equal to er just above it
    er(i_d,1:NPML_y+1+n_os_y)=er(i_d,NPML_y+1+n_os_y);
% top er equal to er just below
    er(i_d,(Ny-NPML_y-(1+n_os_y)):Ny)=er(i_d,Ny-NPML_y-1-n_os_y);
% bottom mr equal to mr just above
    mr(i_d,1:NPML_y+1+n_os_y)=mr(i_d,NPML_y+1+n_os_y);
% top mr equal to mr just below
    mr(i_d,(Ny-NPML_y-(1+n_os_y)):Ny)=mr(i_d,Ny-NPML_y-1-n_os_y);  
    end
    
% define y=mx+b to define a sigma that slowly varies from 0 to
% sigma max, to decay wave
    slope_x=(sigma_1-sigma_PML_max)/(NPML_x-1);
    slope_y=(sigma_1-sigma_PML_max)/(NPML_y-1);
    b_x=sigma_1-slope_x*NPML_x;
    b_y=sigma_1-slope_y*NPML_y;
% set left and right PML sigma values
for i=1:NPML_x    
        sigma_e(i,:)=slope_x*i+b_x;        
        sigma_e(Nx-i+1,:)=slope_x*i+b_x;
        sigma_m(i,:)=m_fac_left.*sigma_e(i,:);
        sigma_m(Nx-i+1,:)=m_fac_right.*sigma_e(Nx-i+1,:);    
end
% set top and bottom PML values
for j=1:NPML_y
        sigma_e(:,j)=slope_y*j+b_y+sigma_e(:,j);        
        sigma_e(:,Ny-j+1)=slope_y*j+b_y+sigma_e(:,Ny-j+1);
        sigma_m(:,j)=m_fac_bottom.*sigma_e(:,j);
        sigma_m(:,Ny-j+1)=m_fac_top.*sigma_e(:,Ny-j+1);    
end    
end

