
clear all
close all
clc

% This code solves for 2D force density on a Bragg grating from a TM mode
% right now solves for each model in serial, should probably change that to
% parallel

addpath TM_EM_functions        % folder with sub-functions
addpath material_data       % folder with material data
addpath position_functions
addpath .\Global_Functions
addpath .\Global_Functions\position_functions
addpath .\Global_Functions\gui_functions

% Make Folder Directories
mkdir('.\geom_data');
mkdir('.\field_data');
mkdir('.\exported_png');

% Define Models
%models=[ 'AB ';'MN ';'AMP';'EL ';'Chu'];

models=[ 'MN '];

%models=[ 'AB ';'AMP';'EL '];


%models=['EL '];



models = cellstr(models);

% Define sim case, where geometry will change
sim_case='grating';             
save_mode=0;                                % Turn on/off saveing the whole work space
write_data=0;                               % Turn on/off writing to text fiels

file_prefix='TM_prelim';

for sim_j=1:length(models)

% close all  
    keep file_prefix models sim_j save_mode sim_case index_n N_run y_disp write_data fx_peak_models fy_peak_models      % clear everything except which model to use.

    clf                     % Clear models inbetween simulations
    plot_on=1;              % turn plotting on or off
    dispersion_on=0;
    model=models(sim_j);    % Define the current model
    model_name=char(model); % store as charecter array for string writing

    % FIGURE INITIALIZATION

% Display figures as
%___________
%| p1 p2 p3|
%| p4 p5 p6|



    s_d=get(0,'ScreenSize');        % Screen [0 0 width height]
    sw=s_d(3);                      % Screen width
    sh=s_d(4);                      % Screen height


    p1=[0 sh/2 sw/3 sh/2];          % Original Figure, top left
    p2=right_of(p1,p1(3),p1(4));    % to the right of p1
    p3=right_of(p2,p1(3),p1(4));    % to the right of p3
    p4=below_of(p1,p1(3),p1(4));    % below p1, second row
    p5=right_of(p4,p1(3),p1(4));    % to the right of p4
    p6=right_of(p5,p1(3),p1(4));    % to ythe right of p5


    fig_1=figure(1);                
    set(fig_1,'name','Source Profile')   
    set(fig_1,'DefaulttextFontSize',14)
    set(fig_1,'OuterPosition',p1)

    fig_2=figure(2);
    set(fig_2,'name','Hz field')
    set(fig_2,'DefaulttextFontSize',14)
    set(fig_2,'OuterPosition',p2)

    fig_3=figure(3);
    set(fig_3,'name','X-Momentum')
    set(fig_3,'DefaulttextFontSize',14)
    set(fig_3,'OuterPosition',p3)

    fig_4=figure(4);
    set(fig_4,'name','Y-momentum')
    set(fig_4,'DefaulttextFontSize',14)
    set(fig_4,'OuterPosition',p4)

    fig_5=figure(5);
    set(fig_5,'name','<Normal force>')
    set(fig_5,'DefaulttextFontSize',14)
    set(fig_5,'OuterPosition',p5)

    fig_6=figure(6);
    set(fig_6,'name','<Lateral force>')
    set(fig_6,'DefaulttextFontSize',14)
    set(fig_6,'OuterPosition',p6)

    

%% SI UNITSfx_avg_
    meters= 1;
    nm=meters*1e-9;
    micrometers=1E-6;
    milimeters=1E-3;
    femptoseconds=1e-15;
    fs=1e-15;
    mu_o=4*pi*10^-7;
    c=299792458;
    eps_o=(1/(c*c*mu_o));
    eta_o=sqrt(mu_o/eps_o);
    Kg=1;


%% Simulation Parameters

    if (strcmp(model,'NA'))
    air_everything=1;               % NA will turn everything to Air
    else
    air_everything=0;
    end
    
%% SIM SETTINGS
    solid_index=0;                  % Turn grating into a solid structure.                  
% Spatial     
    LAMBDA=1064*nm;                 % wavelength of source
%     dx=25*nm;                       % x-step size
%     dy=25*nm;                       % y-step size
    
    n_SiO2=1.434;    %   d=191.39 in 
    n_ZrO2=1.883;    %   d=136.93

%     dx=LAMBDA/10/3;                       % x-step size
%     dy=LAMBDA/10/3;                       % y-step size
%     
    dx=17*nm;                                % x-step size
    dy=dx;                                  % y-step size
    
% Temporal

    dt=(dx/2.0)*1/(c);                  % time step       
    sim_cycles=20;                      % total source cycles to simualte 
    f=c./LAMBDA;                        % Center frequency of pulse
    T=1./f;                             % Center period

    sim_time=(sim_cycles*(1/f));        % total simtime in (s) 
    Nt=round(sim_time/dt);              % number of time steps

    fig_refresh=T;                    % 10 figures per cycle
    
    fig_count=round(fig_refresh/dt);    % time steps per plot
    n_count=0;                          % Plot Counter

    if strcmp(model,'NA')
    end

    
% GEOMETRY

    A_slab=1;                           % Area of cylinder [m^2]
    M_slab=.0001*Kg;                    % Mass of cylinder [Kg]
    M_system=M_slab;                    % Assume mass of system is same slab    
    theta_o=0*2*pi/360;                 % Source Angle
    L_x=2100*nm;                       % x length of object 
    L_y=9500*nm;
    
% If set to 'NA' model, or solid index, the sim does not need to be
% very long in the incident direction 

    if ( strcmp(model,'NA') || solid_index==1)
        L_x=1*LAMBDA;            % x length of object 
        sim_cycles=40;            % total source cycles to simualte
    end
    
% SOURCE GEOMETRY

    gauss_y_width=0.7*L_y;                          % Width of gaussian beam
    gauss_x_width=1*LAMBDA;                         % Width of gaussian beam
    sig_y=1*gauss_y_width/2;                        % Deviation, sigma is half the width                 


    
    spc_x1=0.5*LAMBDA;                              % Space on the left before object
    spc_x2=5*dx;                                    % Space on the right after the object
    spc_y1=.2*LAMBDA;                               % Bottom  (y=0 and grater)              
    spc_y2=.2*LAMBDA ;                              % Space between slab and Y_PML
    
    source_type='cw';                               % source propagation direction 
    source_direction=1;                             % 1 for forwards. -1 for backwards
    top_hat_profile=1;                              % turn on or off top hat profile
    straight_profile=0;
    gauss_t_on=0;                                   % multiply source by a time gaussian, essentially
                                                    % Making it a pulse
    rise_t_on=1;                                    % turn on a rise time for pulse
    
% GRID SIZE

    x_size=spc_x1+spc_x2+L_x;                       % x[m] size of simulation
    y_size=max([L_y gauss_y_width])+spc_y1+spc_y2;  % y[m] size of simulation
    
    NPML=.5*LAMBDA;                                 % PML size
    NPML_x=round(NPML/dx);                          % size of PML along x    
    NPML_y=NPML_x;                                  % size of PML along y    
        
    Nx=round(x_size/dx)+2*NPML_x;                   % Grid x length
    Ny=round(y_size/dy)+2*NPML_y;                   % Grid y length
    
    slab_time=(spc_x1)./(c);                      % time it takes to reach the slab
    %avg_time=1.5*x_size/(c/1.88);                          % start taking average 9 cycles before simulation is finished   
    
    avg_time=(sim_cycles-10)*T;
    
    x=[0:1:Nx-1]*dx;                        % x axis 
    y=[0:1:Ny-1]*dy;                        % y axis 
    t=[0:1:Nt-1]*dt;                        % time axis in [s]   

    c_x=mean(x);                            % center x location of axis
    c_y=mean(y);                            % center y of grid used before rotation?
    c_y_object=mean(y);                     % y location of slab
    c_x_object=NPML_x*dx+spc_x1+L_x/2;      % c_x of object

    gauss_x_avg=c_x;                        % source gauss_x_avg
    gauss_y_avg=c_y;                        % source gauss_y_abg


    
%% DEFINE INDEX OF REFRACTION

    er_background=1;                             % Background \eps_r
    mr_background=1;                             % Background \mu_r      
    er=er_background.*ones(Nx,Ny);
    mr=mr_background.*ones(Nx,Ny); 

% refractive index's
% SiO2   n= 1.434 (1000nm)
% Zn0_2  n= 1.883 (      ) 

    er_SiO2=1.434^2;    %   d=191.39 in 
    er_ZrO2=1.883^2;    %   d=136.93

% Define first layer

    er_1=er_SiO2;                       % Sio2 at 1064=1.45
    mr_1=1;     

    % d_1=LAMBDA/(4*sqrt(er_1));        % formula for ideal quarte wave plate
    d_1=191.39*nm;                      % actual size from manufacturer

    er_coating=er_SiO2;
    d_coating=191.39*nm;

    gamma_e_1=0;                        % electric drude damping f
    omega_e_1=0;                        % electric drude resonance f

    gamma_m_1=0*5E14;                   % magnetic damping frequency
    omega_m_1=0*6E15;                   % magnetic resonance frequency

    sigma_e_1=0;                        % Electrical Conductivity
    sigma_m_1=0;                        % Magnetic Conductivity

% Define second layer
    d_2=136.93*nm;                      % Size of first layer as a quarter of the wavelength
    er_2=er_ZrO2;                       % Permittivity from tomaz ZrO2
    if solid_index==1
      er_2=er_1;  
    end
    mr_2=1;  

% d_2=LAMBDA/(4*sqrt(er_2));            % formula for ideal quarte wave plate

    gamma_e_2=0;                        % Electric drude damping                            
    omega_e_2=0;                        % Magnetic drude resonance

    gamma_m_2=0;                        % Damping frequency of object
    omega_m_2=0;                        % Electric plasma frequency of object

    sigma_e_2=0;                        % Electrical Conductivity
    sigma_m_2=0;                        % Magnetic Conductivity  

% Define Bragg Grating Layers within box
% with corner locations defined by cx1,cx2,cy1,cy2
    [er,cx1,cx2,cy1,cy2,box_1x_mat,box_1y_mat,box_2x_mat,box_2y_mat]=...
        create_bragg_grating_v2(L_x,L_y,c_x_object,c_y,er_1,er_2,d_1,d_2,er,dx,dy); 

    n_coating=round(d_coating/dx)+1; % index of coating start    
    er(cx1-n_coating:cx1,cy1:cy2)=er_coating;

    coating_box_x=[x(cx1-n_coating) x(cx1-n_coating) x(cx1) x(cx1) x(cx1-n_coating)];
    coating_box_y=[y(cy1) y(cy2) y(cy2) y(cy1) y(cy1)];

   source1_x=round((x(NPML_x)+.5*(x(cx1-n_coating)-x(NPML_x)) )/dx);  % x-location of source, if pointed along x
    source1_y=round(Ny/2);                  % y-location of source  if pointed along y


    % Write geometry locations to text files

    
    if write_data==1

    write_dat(strcat('./geom_data/er_coating_box_',sim_case,'.dat'),coating_box_x./nm,coating_box_y./nm);
    write_dat(strcat('./geom_data/er_1_boxes_',sim_case,'.dat'),x(box_1x_mat)./nm,y(box_1y_mat)./nm);
    write_dat(strcat('./geom_data/er_2_boxes_',sim_case,'.dat'),x(box_2x_mat)./nm,y(box_2y_mat)./nm);
    write_dat(strcat('./geom_data/grid_',sim_case,'.dat'),[x(1) x(end) x(end) x(1) x(1)]*1e6,[y(end) y(end) y(1) y(1) y(end)]*1e6);
    write_dat(strcat('./geom_data/pml_',sim_case,'.dat'),[x(NPML_x) x(Nx-NPML_x) x(Nx-NPML_x) x(NPML_x) x(NPML_x)]*1e6,[y(Ny-NPML_y) y(Ny-NPML_y) y(NPML_y) y(NPML_y) y(Ny-NPML_y)]*1e6);


    end

    material_location_matrix=1.*(er==er_1)+1.*(er==er_2);   % Matrix of 1's wherever there is material
    slab_mat=material_location_matrix;                      % Designate "Slab" as the object being pushed

% Define other material coefficients
    gamma_e=gamma_e_1.*(er==er_1)+gamma_e_2.*(er==er_2);
    gamma_m=gamma_m_1.*(er==er_1)+gamma_m_2.*(er==er_2);
    omega_m=omega_m_1.*(er==er_1)+omega_m_2.*(er==er_2);
    omega_e=omega_e_1.*(er==er_1)+omega_e_2.*(er==er_2);

    if ((strcmp(model,'AB'))||(strcmp(model,'MN')))

        [Fx,Fy]=gradient(material_location_matrix);         %  take a gradient
        mat_1=1.*(Fx~=0);                                   % wherever there is a transition, add to integration location
        mat_2=1.*(Fy~=0);                                   % this will capture boundary differences, needed for AB and MN                                   
        slab_mat=material_location_matrix+(mat_1+mat_2);
        slab_mat=1.*(~(slab_mat==0));                       % set to one wherever material exists

% Second gradient, averageing will sometimes look i-2 spatial steps away   

        [Fx,Fy]=gradient(slab_mat);
        mat_1=1.*(Fx~=0);
        mat_2=1.*(Fy~=0);
        slab_mat=slab_mat+mat_1+mat_2;
        slab_mat=1.*(~(slab_mat==0));
    end
    
% set system mat as everything else that is not the slab       
    system_mat=1.*(~(slab_mat==1));

% Define x and y axis for slab_mat
% to call surf(x,y,slab_mat)

    [I_slab J_slab]=find(1.*(slab_mat==0)); % Find locations where slab exists

% find corner nodes
    x_slab_min=min(I_slab)*dx;
    x_slab_max=max(I_slab)*dx;
    y_slab_min=min(J_slab)*dy;
    y_slab_max=max(J_slab)*dy;
% create matrix for plotting
    x_slab=[x_slab_min:dx:x_slab_max];
    y_slab=[y_slab_min:dy:y_slab_max];

 
 % if air_everything is on, 
 % turn everything into 1
 if air_everything==1
   
    er=ones(Nx,Ny);
    mr=er;
    gamme_e=zeros(Nx,Ny);
    gamme_m=zeros(Nx,Ny);

    omega_m=zeros(Nx,Ny);
    omega_e=zeros(Nx,Ny);
    
end
 
 
%% CREATE PML
% Values are found by trial and error
% can use decay depth equation as well
    sigma_1=0;                          % Minimum PML sigma value at boundary of simulation space with PML
    sigma_PML_max=90000;                % Maximum PML sigma value at edge of PML
    n_os_x=round(spc_x2/dx)+1;          % Desired x offset to "match" any permitivity next to the PML
    n_os_y=round(2*spc_y1/dx);          % y offset

% Call function to place PML values of er,mr, and sigma's on the gird
[ er,mr,sigma_e,sigma_m ] = create_PML(sigma_1,sigma_PML_max,er,mr,NPML_x,NPML_y,n_os_x,n_os_y );

%% Source
    E_pulse= 160/1000;                              % mJ of pulse
% average_power=.0016;                              % average power determined from "average_power_calculate.m")
    average_power=6.93E6;                           % found by trail and error?
    P_o=1;                                          % momentum normalization, 1 to report absolute values.
    r_puck=2.058*milimeters;

%power_density=average_power/(pi*r_puck^2);    % W/m^2 of the pulse <S>
power_density=6.4755E11;

    a=sqrt(power_density*(2*eta_o));                % amplitude of electric field for <S>
   
   % a=a*sqrt(1.095); % discrepency between TE and TM 
    
    Sx_theory=dy/dx*a^2/(2*eta_o);                  % theoretical plane wave <S> for comparions

% TIME PROFILE
    t_avg=(4*T);                                    % Location of pulse peak in time     
    sig_t=1*(1/f);                                  % Temporal gaussian width of pulse    
    gauss_t=1*exp(-1.*((t-t_avg).^2)/(2*sig_t^2));  % Time Gaussian  

    if gauss_t_on==0                                % Turn of time dependance if needed
    gauss_t=ones(size(gauss_t));
    end

    if rise_t_on==1  
    tau=5*T;                                                        % Rise time of pulse
    gauss_t=(1-1.*exp(-t./(sig_t))).*gauss_t; 
    end

    gauss_y=gauss_create(gauss_y_avg,sig_y,[0:Ny-1]*dy);            % Define spatial gaussian
    gauss_y=gauss_y./max(gauss_y);

% Define index to the left and right bassed on gaussiean width
% assume this width is where the spatial profile is non-zero
    cy1_profile=round(Ny/2)-round(gauss_y_width/(2*dy));
    cy2_profile=round(Ny/2)+round(gauss_y_width/(2*dy));

    if straight_profile==1
       gauss_y=ones(size(gauss_y)); 
    end
    
    if top_hat_profile==1 

        alpha_c=.2;                                                     % Alpha value of tukey window
        gauss_y=sqrt(tukey_window_create(gauss_y_width,alpha_c,y,dy));  % call function to create tukey window based on pulse_approximation.pdf 
        % interpolate from a picuter of eta(r) to plot FDTD's eta(r)
        % against the one from pulse_appoximation.pdf

        [tomaz_r,tomaz_eta]=picture_to_plot('./Resources/Pulse_Approximation','png',0,2.1,0,1.05,[63 61 153 ],200);
        tomaz_eta=tomaz_eta./max(tomaz_eta);           % normalize if needed
        tomaz_r=[fliplr(-tomaz_r) tomaz_r(2:end)];     % define over
        tomaz_eta=[fliplr(tomaz_eta) tomaz_eta(2:end)]-min(tomaz_eta);
        Ny_gauss=length(cy1_profile:cy2_profile);
        Nr_tomaz=length(tomaz_r);
        tomaz_eta=interp1([0:1:Nr_tomaz-1]./(Nr_tomaz-1),tomaz_eta,[0:Ny_gauss-1]./(Ny_gauss-1));


    end

    if write_data==1
        write_dat('./geom_data/gauss_y_profile_TM.dat',y./nm,gauss_y)
    end


%% Output Fields
% Initialize Fields and Current matricies
    Hz=zeros(Nx,Ny);
    Bz=zeros(Nx,Ny);
    Sx=zeros(Nx,Ny);

    Sx_flux=zeros(1,Nt);            % instantaneous power
    Sx_flux_ref=zeros(1,Nt);        % reflected power
    Sx_flux_ref_avg=zeros(1,Nt);    % average reflected power

    Sx_avg=zeros(1,Nt);             % average incident power

 % Spatial average quantities
 
    Hz_at_x=zeros(Nx,Ny);           % Hz placed at Ex
    Hz_at_y=zeros(Nx,Ny);           % Hz placed at Ey
    Bz_at_x=zeros(Nx,Ny);           % Bz placed at Ex
    Bz_at_y=zeros(Nx,Ny);           % Bc placed at Ey
    
% Initialize field values

    Ex=zeros(Nx,Ny);
    Ey=zeros(Nx,Ny);
    Dx=zeros(Nx,Ny);
    Dy=zeros(Nx,Ny);

% Dispersive current densities

    Jxd=zeros(Nx,Ny);
    Jyd=zeros(Nx,Ny);
    Jmzd=zeros(Nx,Ny);
        
% Static current densities

    Jxs=zeros(Nx,Ny);
    Jys=zeros(Nx,Ny);
    Jmzs=zeros(Nx,Ny);
    
% Dispersive polarization and magnetization.

    Pxd=zeros(Nx,Ny);
    Pyd=zeros(Nx,Ny);
    Mzd=zeros(Nx,Ny);
    
% Static Polarization and magnetizaton
    
    Pxs=zeros(Nx,Ny);
    Pys=zeros(Nx,Ny);
    Mzs=zeros(Nx,Ny);
    
% Total polarization and currents
    
    Px=zeros(Nx,Ny);
    Py=zeros(Nx,Ny);
    Mz=zeros(Nx,Ny);
    Jmz=zeros(Nx,Ny);
    Jex=zeros(Nx,Ny);
    Jey=zeros(Nx,Ny);


% Enclosure/"system" variables
    system_momentum_x=zeros(1,Nt);
    system_momentum_y=zeros(1,Nt);
% COM
    x_bar_system=zeros(1,Nt);
    y_bar_system=zeros(1,Nt);
% Force
    F_system_x=zeros(1,Nt);
    F_system_y=zeros(1,Nt); 
                
    system_acceleration_x=zeros(1,Nt);    %[m/s^2]        Total acceleration
    system_velocity_x=zeros(1,Nt);        %[m/s]          Total Velocity        
    system_displacement_x=zeros(1,Nt);    %[m]            Displacement

    system_acceleration_y=zeros(1,Nt);    %[m/s^2]        Total acceleration      
    system_velocity_y=zeros(1,Nt);        %[m/s]          Total Velocity        
    system_displacement_y=zeros(1,Nt);    %[m]            Displacement

    x_bar_system_o=mean(x);
    x_bar_system_contribution=zeros(1,Nt);
    y_bar_system_o=mean(y);
    y_bar_system_contribution=zeros(1,Nt);
   
% Total variables
    Tx=zeros(Nx,Ny);                                % Tx component
    Ty=zeros(Nx,Ny);                                % Ty component
    G_x=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_x_n_prev=0;
    
    g_mech_x=zeros(Nx,Ny);                          % Momentum density
    g_mech_y=zeros(Nx,Ny);                          % Momentum density

    x_bar_total=zeros(1,Nt);    
    y_bar_total=zeros(1,Nt);  
                
% Pulse_Output_Variables
    W=zeros(Nx,Ny);                     %[J/m^3]        Energy Density
    M_pulse=zeros(1,Nt);                %[Kg]           Total pulse mass
    pulse_energy=zeros(1,Nt);           %[J]            Total Energy
    
    source_spatial_profile=zeros(1,Ny);

    % pulse x-variables
    pulse_momentum_x=zeros(1,Nt);               %[kg(m/s)]      Total Momentum


    x_bar_pulse=zeros(1,Nt);            %[m]            x- center of mass
    x_bar_pulse_contribution=zeros(1,Nt);%[m]per kg      % contribution of ceneter of mass

    % y-variables
    pulse_momentum_y=zeros(1,Nt);       %[kg(m/s)]      Total Momentum

    G_y=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_y_n_prev=0;
    y_bar_pulse=zeros(1,Nt);            %[m]            x- center of mass
    y_bar_pulse_contribution=zeros(1,Nt);           %[m]per kg      % contribution of ceneter of mass

       
% Slab Output variables


% x-variables
    F_slab_x=zeros(1,Nt);
    fx_slab_middle=zeros(Nx,1);
    fx_slab_middle_avg=zeros(Nx,1); % cumulative average 
    f_n_count=1;                    % counter to start cumulative average
    f_n_count_2d=1;                 % counter for 2d cumulative average
    F_slab_y=zeros(1,Nt);
    F_slab_x_avg=zeros(1,Nt);
    F_slab_y_avg=zeros(1,Nt);
    fy_slab=zeros(size(slab_mat));
    fx_slab=zeros(size(slab_mat));

    fy_slab_avg=zeros(size(slab_mat));
    fx_slab_avg=zeros(size(slab_mat)); % avg over x-y profile


    fx_avg_sum_track=zeros(1,Nt);
    fy_avg_sum_track=zeros(1,Nt);
    

    % Track 2d force over 10 cycles
%     moving_average_window_size=round(10*T/dt);
%     fx_slab_moving_average=zeros(Nx,Ny, moving_average_window_size);
%     fy_slab_moving_average=zeros(Nx,Ny, moving_average_window_size);


% intialize int
    fy_int_avg=zeros(1,length(y_slab));
    fx_int_avg=zeros(1,length(y_slab));

    slab_momentum_x=zeros(1,Nt);        % [kg(m/s)]      Total Momentum  
    diff_momentum_x=zeros(1,Nt);

    slab_acceleration_x=zeros(1,Nt);    % [m/s^2]        Total acceleration        
    slab_velocity_x=zeros(1,Nt);        % [m/s]          Total Velocity        
    slab_displacement_x=zeros(1,Nt);    % [m]            Displacement

    x_bar_slab_o=(1/(sum(sum(slab_mat))*dx*dy)).*sum(x(1:Nx).*sum((slab_mat(1:Nx,1:Ny))'.*dy)).*dx;

    x_bar_slab=zeros(1,Nt);             % [m]            x-center of mass
    x_bar_slab_contribution=zeros(1,Nt);% [m]per kg      % contribution of ceneter of mas


% y-variables

        
    slab_momentum_y=zeros(1,Nt);            %[kg(m/s)]      Total Momentum    

    slab_acceleration_y=zeros(1,Nt);        %[m/s^2]        Total acceleration
    slab_velocity_y=zeros(1,Nt);            %[m/s]          Total Velocity
    slab_displacement_y=zeros(1,Nt);        %[m]            Displacement
    y_bar_slab=zeros(1,Nt);                 %[m]            x-center of mass
    y_bar_slab_contribution=zeros(1,Nt);    %[m]per kg      % contribution of ceneter of mas

    y_bar_slab_o=(1./(sum(sum(slab_mat))*dx*dy)).*sum(y(1:Ny).*sum((slab_mat(1:Nx,1:Ny)).*dy)).*dx;

 %% Extra output variables
 
    total_momentum_x=zeros(1,Nt);
    total_momentum_y=zeros(1,Nt);
    gx_center_prev=zeros(1,length(source1_x:cx2+4));
    gx_center=zeros(1,length(source1_x:cx2+4));
    fx_center=zeros(length(gx_center),Nt);


%% COEFFICIENTS

    
    H1=(-mu_o.*mr./dt+sigma_m./2)./(-mu_o.*mr./dt-sigma_m./2);        % PML only
    H2=1./(-mu_o.*mr./dt-sigma_m./2);
    
      E1=(eps_o.*er./dt-sigma_e./2)./(eps_o.*er./dt+sigma_e./2);        % PML only
    E2=1./(eps_o.*er./dt+sigma_e./2)/dx;

    

%% Indexes

    i=4:Nx-4;
    j=4:Ny-4;
    
    % index to take average att
    
    [val, peak_index]=max(diff(gauss_y));
    
    %peak_index=round(Ny/2)-round(LAMBDA/4/dy);
    
    % peak index should be 420 next time to get proper reading...
% moving_average_counter=1;
% plot_moving_average=0;
% index locations to track incident source power
    srci=round(.5*(x(source1_x)+x(cx1))/dx);                            % Track halfway between source and first interface cx1
    srcj=round(Ny/2)-round(L_y/(4*dy)):round(Ny/2)+round(L_y/(4*dy));   % only capture 1/2 L_y centered ?

% index location to track reflected source
    refi=source1_x-2;                               % look just behind the source                             
    refj=srcj;

% Figure 1, Source Profile
if top_hat_profile==1
    if plot_on==1
    set(fig_1,'name','Source Profile')
    figure(1)
    h_source_input=plot(y,gauss_y,'color','r');
    hold on
    h_tomaz=plot(y(cy1_profile:cy2_profile),tomaz_eta,'--black');
    hold on
    % h_cy1=line([y(cy1) y(cy1)],[0 1],'color','black');
    % h_cy2=line([y(cy2) y(cy2)],[0 1],'color','black');
    h_NPML_1=line([y(NPML_y) y(NPML_y)],[0 1],'color','b','linestyle','--');
    h_NPML_2=line([y(Ny-NPML_y) y(Ny-NPML_y)],[0 1],'color','b','linestyle','--');
    hold on
    h_source_fdtd=plot(y,source_spatial_profile,'--g');
    legend([h_tomaz,h_NPML_1,h_source_input,h_source_fdtd],'Laser Pulse.pdf','NPML_y','\eta{(r)}','fdtd \eta{(r)}','location','south')

    end
    end
% PLOT HZ FIELD
    figure(2)                           
    h_Hz=surf(x*1e6,y*1e6,Hz');                               % surf the Hz field
    hold on
    [dd, h_er]=contour(x*1e6,y*1e6,er',3,'color','black');               % Plot geometry
    hold on
    % Draw PML limes
    line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %  
    hold on
    % Line where source is measured
    h_src_meas=line(srci*1e6*[dx dx],[srcj(1) srcj(end)]*dy*1e6,[1.1 1.1],'color','m');
    hold on
    % line where reflection is measured
    h_ref=line(refi*1e6*[dx dx],[refj(1) refj(end)]*dy*1e6,[1.1 1.1],'color','cyan');
    % line where souce is produced
    h_source=line(source1_x*1e6*[dx dx],[refj(1) refj(end)]*dy*1e6,[1.1 1.1],'color','green');

    xlabel('x axis [{\mu}m]') % 
    ylabel('y axis [{\mu}m]') % 
    shading flat
    %caxis([-2 2]*cmax)      
    legend([h_src_meas,h_ref,h_source],'src meas.','ref meas.','source','location','southeast');       
    drawnow
    view([0 90])  
    title(model)
% FIGURE 2 X-MOMENTUM
set(0, 'CurrentFigure', fig_3)
    h_Px_pulse=plot(t,pulse_momentum_x,'color','r');
    hold on
    h_Px_slab=plot(t,slab_momentum_x,'color','b');
    hold on
    h_Px_system=plot(t,system_momentum_x,'--','color','black');
    hold on
    h_Px_total=plot(t,(slab_momentum_x+pulse_momentum_x+system_momentum_x),'color','g');
    hold on
    title(model)
    h_avg_line=plot(0,0);
    xlabel('time(s)')
    ylabel(' P_x [kgm/s]')
    legend('Pulse','Object','System','Total','Location','northwest')
    
% FIGURE 4 Y-MOMENTUM
set(0, 'CurrentFigure', fig_4)
    h_Py_pulse=plot(t,pulse_momentum_y,'color','r');
    hold on
    h_Py_slab=plot(t,slab_momentum_y,'color','b');
    hold on
    h_Py_system=plot(t,system_momentum_y,'--','color','black');
    hold on
    h_Py_total=plot(t,(slab_momentum_y+pulse_momentum_y+system_momentum_y),'color','g');

    title(model)
    xlabel('time(s)')
    ylabel(' P_y [kgm/s]')
    legend('Pulse','Object','System','Total','Location','northwest')
% Figure 5 force profiles
set(0, 'CurrentFigure', fig_5)
subplot(2,1,1)
    h_fx_surf=surf(x_slab*1e6,y_slab*1e6,fx_slab_avg');
        hold on
    [dd, h_er_fx]=contour(x*1e6,y*1e6,er',3,'color','black');               % Plot geometry
    hold on
    % Draw PML limes
    line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %  
    view([0 90])
    shading flat
subplot(2,1,2)
    h_fx_avg=plot(0,0,'blue');
    hold on
        h_fx_avg_moving=plot(0,0,'r--');

    xlabel('x')
    ylabel('<fx>')
    legend('sum','moving window')

% Figure 6 lateral force profile 
set(0, 'CurrentFigure', fig_6)
subplot(2,1,1)
    h_fy_surf=surf(x_slab*1e6,y_slab*1e6,fy_slab_avg');
    hold on
    [dd, h_er_fy]=contour(x*1e6,y*1e6,er',3,'color','black');               % Plot geometry
    hold on
    % Draw PML limes
    line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %  
    view([0 90])
    shading flat
subplot(2,1,2)
    h_fy_avg=plot(0,0,'blue');
    hold on
    h_fy_avg_moving=plot(0,0,'r--');
    view([0 90])
    shading flat
     ylabel('<fy>')

         legend('sum','moving window')

% clear uneccessary values for the rest of the code
clear gamma_e gamma_m  sigma_m sigma_e A1 A2 X1 X2 E omega_e omega_m
clear FWHM E_source H_source E_y E_tot 
clear N NPML N_y_gauss W_source b_x b_y den_x den_y
clear eps_1 eps_2  er_Im_Si er_Re_Si  
clear f_end f_start gauss_y_end gauss_y_start gauss_x_avg 
clear gamma_e_material gamma_m_material mu_2 n_Si_exp
clear k_Si_exp k_mat lam0 lambda_Si norm_G omega_e 
clear m_fac omega_m 
clear q s_d sh sig_t sig_y sigma_PML_e sigma_PML_max sim_cycles sim_time
clear slope_x slope_y  sw t_avg
clear mr_d  er_d
clear mat_1 mat_2
    for n= 1:Nt
        
% Store Previous values 
% Previous Ex field values
    Ex_n_prev=Ex;
    Dx_n_prev=Dx;   
    Jxd_n_prev=Jxd(2:Nx,2:Ny);      
    Jxs_n_prev=Jxs(2:Nx,2:Ny);      
% Previsous Ey field
    Ey_n_prev=Ey;
    Dy_n_prev=Dy;
    Jyd_n_prev=Jyd(2:end,2:end);   
    Jys_n_prev=Jys(2:end,2:end);   
% Previous Hz field values
    Hz_n_prev=Hz;
    Bz_n_prev=Bz;
    Jmz_n_prev_fr=Jmz;
    Jmzd_n_prev=Jmzd;    
% Previous force values
    g_mech_x_n_prev=g_mech_x;
    g_mech_y_n_prev=g_mech_y;
   


%% Update Ex,Ey (n+1)
 % simple Ex 
% Ex(2:Nx,2:Ny)=Ex(2:Nx,2:Ny)+E1(2:Nx,2:Ny).*(Hz(2:Nx,2:Ny)-Hz(2:Nx,1:Ny-1)); % lossless


% Update Ex
if dispersion_on==1
    Ex(2:Nx,2:Ny)=(X3).*Ex(2:Nx,2:Ny)-X4.*Jxd(2:Nx,2:Ny)+X5.*(Hz(2:Nx,2:Ny)-Hz(2:Nx,1:Ny-1));
    
    Jxd(2:Nx,2:Ny)=P1.*Jxd(2:Nx,2:Ny)+P2.*(Ex(2:Nx,2:Ny)+Ex_n_prev);        % Dispersive current density
    Pxd(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+(dt/2).*(Jxd(2:Nx,2:Ny)+Jxd_n_prev);

    Jxs(2:Nx,2:Ny)=-Jxs(2:Nx,2:Ny)+(2/dt)*C1.*(Ex(2:Nx,2:Ny)-Ex_n_prev);    % Static current density
    Pxs(2:Nx,2:Ny)=Pxs(2:Nx,2:Ny)+(dt/2)*(Jxs(2:Nx,2:Ny)+Jxs_n_prev);

    Jex(2:Nx,2:Ny)=Jxs(2:Nx,2:Ny)+Jxd(2:Nx,2:Ny);                           % Total current density
    Px(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+Pxs(2:Nx,2:Ny);
    Dx=eps_o.*Ex+Px;
     

else
   Ex(2:Nx,2:Ny)=E1(2:Nx,2:Ny).*Ex(2:Nx,2:Ny)+E2(2:Nx,2:Ny).*(Hz(2:Nx,2:Ny)-Hz(2:Nx,1:Ny-1)); %PML only
 
  Dx=eps_o.*er.*Ex;
end

% UPDATE Ey
% simplest form % Ey(2:Nx,2:Ny)=Ey(2:Nx,2:Ny)-E1(2:Nx,2:Ny).*(Hz(2:Nx,2:Ny)-Hz(1:Nx-1,2:Ny));

if dispersion_on==1
    Ey(2:Nx,2:Ny)=Y1.*Ey(2:Nx,2:Ny)-Y2.*Jyd(2:Nx,2:Ny)-Y3.*(Hz(2:Nx,2:Ny)-Hz(1:Nx-1,2:Ny));

    Jyd(2:Nx,2:Ny)=P1Y.*Jyd(2:Nx,2:Ny)+P2Y.*(Ey(2:Nx,2:Ny)+Ey_n_prev);
    Pyd(2:end,2:end)=Pyd(2:end,2:end)+(dt/2)*(Jyd(2:end,2:end)+Jyd_n_prev);

    Jys(2:Nx,2:Ny)=-Jys(2:Nx,2:Ny)+(2/dt)*C1.*(Ey(2:Nx,2:Ny)-Ey_n_prev);
    Pys(2:end,2:end)=Pys(2:end,2:end)+(dt/2)*(Jys(2:end,2:end)+Jys_n_prev);

    Jey(2:Nx,2:Ny)=Jys(2:Nx,2:Ny)+Jyd(2:Nx,2:Ny);
    Py(2:end,2:end)=Pyd(2:end,2:end)+Pys(2:end,2:end);
    Dy=eps_o.*Ey+Py;

else
 Ey(2:Nx,2:Ny)=E1(2:Nx,2:Ny).*Ey(2:Nx,2:Ny)-E2(2:Nx,2:Ny).*(Hz(2:Nx,2:Ny)-Hz(1:Nx-1,2:Ny));
   Dy=eps_o.*er.*Ey;  
   
end
%% FDTD Hz Source
for gi=1:length(gauss_y)  
Hz(source1_x,gi)=gauss_t(n).*gauss_y(gi).*(a/(2*eta_o)).*sin(2*pi*f*(t(n)))+Hz(source1_x,gi);  % X- Guassian %         
end             


%% W,G,S,T Calculate FIVE FORMS

if (strcmp(model,'MN'))
%[W]=calculate_W_MN(c^2,i_W,j_W,Dx,Dx_n_prev,Dy,Dy_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);
% turn off W calculation for speed
    W(source1_x-1:source1_x+1,gi)=0;                                   % Remove energy at source location

    Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));
  Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
  
  % Bring Bz in space to Dx, but not in time
  
    Bz_at_x(i,j)=(1/2)*(Bz(i,j)+Bz(i,j-1));

  
    G_x_n_prev=G_x;
    G_x(i,j)=Dy(i,j).*Bz_at_y(i,j);                  % G_MN=DXB
    Sx(i,j)=(Dy(i,j).*Bz_at_y(i,j)).*c^2;
    
    G_y_n_prev=G_y;
    % Ty is at Ex
	 
	 % bring Dx to B in time
   
        %Dx_av=(0.5).*(Dx(i,j)+Dx(i,j+1)); % only to place Ty at Bz
        Dx_av=(0.5).*(Dx(i,j)+Dx_n_prev(i,j)); % only to place Ty at Bz

	% Bring Dx to B temporally
   % Gy(i,j)=Ez_at_Bx(i,j).*Hx(i,j)./c^2; TE CODE

   % G_y(i,j)=-1*Dx_av.*Bz_at_x(i,j);   % Gy is at i+1/2,j    
   
   
        G_y(i,j)=-1*Dx(i,j).*Bz_at_x(i,j);   % Ty is at Ex i think    

    % will put tY
[Tx ] = Calculate_Tx_MN( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
   % [Ty ] = Calculate_Ty_MN( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
[Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
		Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );

    
end


if (strcmp(model,'Chu'))    
    %[W]=calculate_W_AB_v2(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);  
     W(source1_x-1:source1_x+1,gi)=0;

    Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
    Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));    

    G_x_n_prev=G_x;        
    G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB  
    Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));

    [ Tx] = Calculate_Tx_Chu(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   

    G_y_n_prev=G_y;
    G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));               
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_Chu( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
end


if (strcmp(model,'EL'))    
      %  [W]=calculate_W_AB_v2(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);
        W(source1_x-1:source1_x+1,gi)=0;

        Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
        Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB    
        [Tx] = Calculate_Tx_EL(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
        Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));
        
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
        % Ty is derivatives about x location 
        % Ty is at i+1/2, and in time at B
      %  [Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL_v3( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
   % Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
    
   [Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL( i,j,Ex,Ey,Dx,...
    Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
end

if (strcmp(model,'AB'))
    
   % [W]=calculate_W_AB_v2(i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W); 
     W(source1_x-1:source1_x+1,gi)=0;

    Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
    Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
    G_x_n_prev=G_x;        
    G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB    
    [ Tx] = Calculate_Tx_AB(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
    Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));

    G_y_n_prev=G_y;
    G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
  
    
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AB( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
		Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
    
end

if (strcmp(model,'AMP'))
    
       % [W]=calculate_W_MN(eps_o*c^2,i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);  
         W(source1_x-1:source1_x+1,gi)=0;

        Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));   
        Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j));
        [Tx] = Calculate_Tx_AMP(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*(Ex(i,j).*Bz_at_x(i,j));                 % 
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AMP( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
       
        
        Sx(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j)).*c^2;
        %Sx(i,j)=mu_o.*(Dy(i,j).*Hz_at_y(i,j)).*c^2;
        G_x(i,j)=Sx(i,j)/c^2;    

end

if (strcmp(model,'NA'))
    
    
       %[W]=calculate_W_MN(eps_o*c^2,i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Bz,Bz_n_prev,dx,dy,dt,W); 
        W(source1_x-1:source1_x+1,gi)=0;

        Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));   
        Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j));
        [Tx] = Calculate_Tx_AMP(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*(Ex(i,j).*Bz_at_x(i,j));                 % 
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AMP( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
        Sx(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j)).*c^2;
        %Sx(i,j)=mu_o.*(Dy(i,j).*Hz_at_y(i,j)).*c^2;
        G_x(i,j)=Sx(i,j)/c^2;    

end


    Sx_flux_ref(1,n)=sum(Sx(refi,refj))./(length(refj));
  Sx_flux(1,n)=sum(Sx(srci,srcj))*dx;
    % If Sx_flux has a value, and more than two periods as passed
    if ((abs(Sx_flux(1,n))>0 )&&(round(t(n)/T)>2))
        
    ni_avg=n-round(T/dt):n;
    Sx_avg(1,n)=mean(Sx_flux(ni_avg)); 
    ref_Sx(1,n)=100*(1-Sx_avg(1,n)./Sx_theory);

    end   
     




%% Update Hz (n+1/2)

    E_term=-1*(1/dx).*(Ey(2:Nx,1:Ny-1)-Ey(1:Nx-1,1:Ny-1))...
            +(1/dy).*(Ex(1:Nx-1,2:Ny)-Ex(1:Nx-1,1:Ny-1));
       
% Hz(1:Nx-1,1:Ny-1)=Hz(1:Nx-1,1:Ny-1)+H1(1:Nx-1,1:Ny-1).*E_term; % lossless

%         
       % simple method
        
if dispersion_on==1
    Hz(1:Nx-1,1:Ny-1)=  ((A3).*Hz(1:Nx-1,1:Ny-1)...    
                        -(A4).*((M3).*Jmzd(1:Nx-1,1:Ny-1))...
                        +(A5).*E_term);

    Jmzd(1:Nx-1,1:Ny-1)=M1.*Jmzd(1:Nx-1,1:Ny-1)+M2.*(Hz(1:Nx-1,1:Ny-1)...
                        +Hz_n_prev(1:Nx-1,1:Ny-1));            

    Jmzs(1:Nx-1,1:Ny-1)=-1*Jmzs(1:Nx-1,1:Ny-1)+(2/dt).*C2(1:Nx-1,1:Ny-1).*(Hz(1:Nx-1,1:Ny-1)-Hz_n_prev(1:Nx-1,1:Ny-1));
    Jmz(1:Nx-1,1:Ny-1)=Jmzs(1:Nx-1,1:Ny-1)+Jmzd(1:Nx-1,1:Ny-1);
        
    Mzd(1:Nx-1,1:Ny-1)=Mzd(1:Nx-1,1:Ny-1)+(dt/2).*(Jmzd(1:Nx-1,1:Ny-1)+Jmzd_n_prev(1:Nx-1,1:Ny-1));
    Mzs(1:Nx-1,1:Ny-1)=C2(1:Nx-1,1:Ny-1).*Hz(1:Nx-1,1:Ny-1);
    Mz=Mzd+Mzs;
   Bz=mu_o.*Hz+Mz;

else 
  
 Hz(1:Nx-1,1:Ny-1)=H1(1:Nx-1,1:Ny-1).*Hz(1:Nx-1,1:Ny-1)-H2(1:Nx-1,1:Ny-1).*E_term; % lossless
% 
% 
     Bz=mu_o.*mr.*Hz;  
end
    
 
for gi=1:length(gauss_y)  
    Ey(source1_x,gi)=source_direction*gauss_t(n).*1*gauss_y(gi).*(.5*a).*sin(2*pi*f*(t(n)))+Ey(source1_x,gi); % Y- Guassian     
end


%% Pulse Calculations
    M_pulse(1,n)=sum(sum(W./c^2))*dx*dy;
    pulse_energy(1,n)=sum(sum(W))*dx*dy;   

    pulse_momentum_x(1,n)=sum(sum(G_x(i,j)))*dx*dy;
    

    
    pulse_momentum_y(1,n)=sum(sum(G_y(i,j)))*dx*dy;

    x_bar_pulse(1,n)=(1/M_pulse(1,n)).*sum(x(i).*sum((W(i,j)./c^2)'.*dy)).*dx;
    y_bar_pulse(1,n)=(1/M_pulse(1,n)).*sum(y(j).*sum((W(i,j)./c^2).*dx)).*dy;
    
%% Slab Calculations

g_mech_x(i,j)=g_mech_x(i,j)-dt.*Tx(i,j)-1*(G_x(i,j)-G_x_n_prev(i,j));
g_mech_y(i,j)=g_mech_y(i,j)-dt.*Ty(i,j)-1*(G_y(i,j)-G_y_n_prev(i,j));    



% stor force density profle along the center
fy_slab=(g_mech_y.*slab_mat-g_mech_y_n_prev.*slab_mat)/dt;
fx_slab=(g_mech_x.*slab_mat-g_mech_x_n_prev.*slab_mat)/dt;

% store cumulative_average
if (n>round(avg_time/dt)) 
    
fy_slab_avg=fy_slab_avg+(fy_slab-fy_slab_avg)./(f_n_count_2d);
fx_slab_avg=fx_slab_avg+(fx_slab-fx_slab_avg)./(f_n_count_2d);
% 
% fx_slab_moving_average(:,:,moving_average_counter)=fx_slab;
%     fy_slab_moving_average(:,:,moving_average_counter)=fy_slab;

%     if moving_average_counter>moving_average_window_size
%         moving_average_counter=1;    
%         fx_slab_moving_average_line=mean(fx_slab_moving_average(source1_x+1:end,round(Ny/2),:),3);
%         fy_slab_moving_average_line=mean(fy_slab_moving_average(source1_x+1:end,peak_index,:),3);
% plot_moving_average=1;
%     end

    

%     moving_average_counter=moving_average_counter+1;


% fy_int_avg=fy_int_avg+(sum(fy_slab_avg,1)*dx-fy_int_avg)./f_n_count_2d;
% fx_int_avg=fx_int_avg+(sum(fx_slab_avg,1)*dx-fx_int_avg)./f_n_count_2d;

fy_int_avg=sum(fy_slab_avg,1)*dx;
fx_int_avg=sum(fx_slab_avg,1)*dx;

fx_avg_sum_track(1,n)=sum(sum(fx_slab_avg))*dx*dy;
fy_avg_sum_track(1,n)=sum(sum(fy_slab_avg))*dx*dy;
%write_dat('./fx_peak.dat',t,fx_avg_peak_track)
%write_dat('./fx_peak.dat',t,fx_avg_peak_track)

f_n_count_2d=f_n_count_2d+1;
% store integration up to surface of materials

end

slab_momentum_x(1,n)=sum(sum(g_mech_x(i,j).*slab_mat(i,j)))*dx*dy;

system_momentum_x(1,n)=sum(sum(g_mech_x(i,j).*system_mat(i,j)))*dx*dy;
    
slab_momentum_y(1,n)=sum(sum(g_mech_y(i,j).*slab_mat(i,j)))*dx*dy;
system_momentum_y(1,n)=sum(sum(g_mech_y(i,j).*system_mat(i,j)))*dx*dy;   
  
if n>1
  
F_slab_x(1,n)=(slab_momentum_x(1,n)-slab_momentum_x(1,n-1))/dt;
slab_acceleration_x(1,n)=1/M_slab*F_slab_x(1,n);     
slab_velocity_x(1,n)=sum(slab_acceleration_x(1,1:n))*dt;
slab_displacement_x(1,n)=sum(slab_velocity_x(1,1:n))*dt;
x_bar_slab(1,n)=x_bar_slab_o+slab_displacement_x(1,n);
 
F_slab_y(1,n)=(slab_momentum_y(1,n)-slab_momentum_y(1,n-1))/dt;
slab_acceleration_y(1,n)=1/M_slab*F_slab_y(1,n);
slab_velocity_y(1,n)=sum(slab_acceleration_y(1,1:n))*dt;
slab_displacement_y(1,n)=sum(slab_velocity_y(1,1:n))*dt;

y_bar_slab(1,n)=y_bar_slab_o+slab_displacement_y(1,n);

F_system_x(1,n)=(system_momentum_x(1,n)-system_momentum_x(1,n-1))/dt;
system_acceleration_x(1,n)=1/M_system*F_system_x(1,n);     
system_velocity_x(1,n)=sum(system_acceleration_x(1,1:n))*dt;
system_displacement_x(1,n)=sum(system_velocity_x(1,1:n))*dt;
x_bar_system(1,n)=x_bar_system_o+system_displacement_x(1,n);
 
F_system_y(1,n)=(system_momentum_y(1,n)-system_momentum_y(1,n-1))/dt;
system_acceleration_y(1,n)=1/M_system*F_system_y(1,n);
system_velocity_y(1,n)=sum(system_acceleration_y(1,1:n))*dt;
system_displacement_y(1,n)=sum(system_velocity_y(1,1:n))*dt;

y_bar_system(1,n)=y_bar_system_o+system_displacement_y(1,n);

% running 2D average
F_slab_x_avg(n)=F_slab_x_avg(n-1)+ (F_slab_x(n)-F_slab_x_avg(n-1) )/n;
F_slab_y_avg(n)=F_slab_y_avg(n-1)+ (F_slab_y(n)-F_slab_y_avg(n-1) )/n;

    
    
    end

	
% y-components    

    
%% System Calculations
M_total=M_slab+M_pulse(1,n)+M_system;
x_bar_total(1,n)= (1/M_total)*(M_slab.*x_bar_slab(1,n)+M_pulse(1,n).*x_bar_pulse(1,n)+M_system*x_bar_system(1,n));
x_bar_pulse_contribution(1,n)=((1/M_total)*(M_pulse(1,n).*x_bar_pulse(1,n)));
x_bar_slab_contribution(1,n)=(1/M_total)*(M_slab.*x_bar_slab(1,n));
 x_bar_system_contribution(1,n)=(1/M_total)*(M_system.*x_bar_system(1,n));
% Y Center of Mass

y_bar_total(1,n)= (1/M_total)*(M_slab.*y_bar_slab(1,n)+M_pulse(1,n).*y_bar_pulse(1,n)+M_system*y_bar_system(1,n));
y_bar_pulse_contribution(1,n)=(1/M_total)*(M_pulse(1,n).*y_bar_pulse(1,n));
y_bar_slab_contribution(1,n)=(1/M_total)*(M_slab.*y_bar_slab(1,n));
y_bar_system_contribution(1,n)=(1/M_total)*(M_system.*y_bar_system(1,n)); 
 

%% Plots
% total variables
total_momentum_x(1:n)=(slab_momentum_x(1:n)+pulse_momentum_x(1:n)+system_momentum_x(1:n));
total_momentum_y(1:n)=(slab_momentum_x(1:n)+pulse_momentum_x(1:n)+system_momentum_x(1:n));




M_tot=M_pulse+M_slab+M_system;
xb_pulse_displacement(6:n)=(M_pulse(6:n)./M_tot(6:n)).*(x_bar_pulse(6:n)-x_bar_pulse(1,6));
xb_slab_displacement(6:n)=(M_slab./M_tot(6:n)).*(x_bar_slab(6:n)-x_bar_slab(1,6));
xb_system_displacement(6:n)=(M_system./M_tot(6:n)).*(x_bar_system(6:n)-x_bar_system(1,6));
xb_total_displacement=xb_pulse_displacement...
                            +xb_slab_displacement...
                            +xb_system_displacement;
   
   
    yb_pulse_displacement(6:n)=(M_pulse(6:n)./M_tot(6:n)).*(y_bar_pulse(6:n)-y_bar_pulse(1,6));
yb_slab_displacement(6:n)=(M_slab./M_tot(6:n)).*(y_bar_slab(6:n)-y_bar_slab(1,6));
yb_system_displacement(6:n)=(M_system./M_tot(6:n)).*(y_bar_system(6:n)-y_bar_system(1,6));
yb_total_displacement=yb_pulse_displacement...
                            +yb_slab_displacement...
                            +yb_system_displacement;


n_count
if n_count==fig_count
  
            if (plot_on==1)

            set(h_Hz,'ZDATA',Hz')

            % UPDATE X-MOMENTUM
            set(h_Px_pulse,'XDATA',t(1:n),'YDATA',pulse_momentum_x(1:n)./P_o);
            set(h_Px_slab,'XDATA',t(1:n),'YDATA',slab_momentum_x(1:n)./P_o);
            set(h_Px_system,'XDATA',t(1:n),'YDATA',system_momentum_x(1:n)./P_o);
            set(h_Px_total,'XDATA',t(1:n),'YDATA',total_momentum_x(1:n)./P_o);
            
            % plot average line
            set(h_avg_line,'XDATA',[avg_time avg_time],'YDATA',[min(slab_momentum_x(1:n)./P_o)...
                max(slab_momentum_x(1:n)./P_o)]);
            % UPDATE Y-MOMENTUM
            set(h_Py_pulse,'XDATA',t(1:n),'YDATA',pulse_momentum_y(1:n)./P_o);
            set(h_Py_slab,'XDATA',t(1:n),'YDATA',slab_momentum_y(1:n)./P_o);
            set(h_Py_system,'XDATA',t(1:n),'YDATA',system_momentum_y(1:n)./P_o);
            set(h_Py_total,'XDATA',t(1:n),'YDATA',total_momentum_y(1:n)./P_o);


            drawnow
            % UPDATE SURFACE PLOT OF FORCE DENSITIES

                    if n>round(avg_time/dt)  
                    set(h_fx_surf,'ZDATA',fx_slab_avg')
                    set(h_fy_surf,'ZDATA',fy_slab_avg')

                    set( h_fx_avg,'XDATA',x(source1_x+2:end),'YDATA',fx_slab_avg(source1_x+2:end,round(Ny/2)))
                    set( h_fy_avg,'XDATA',x(source1_x+2:end),'YDATA',fy_slab_avg(source1_x+2:end,peak_index))

%                             if plot_moving_average==1
%                             set( h_fx_avg_moving,'XDATA',x(source1_x+1:end),'YDATA',fx_slab_moving_average_line)
%                             set( h_fy_avg_moving,'XDATA',x(source1_x+1:end),'YDATA',fy_slab_moving_average_line)
%                             end


                    % Reset figure counter


                    end


            end


n_count=0;
    end
   n_count=n_count+1;
 
    % Data to export at end of simulation 
     


    end

       if write_data==1
         mkdir('./QA_data')
         mkdir('./Figure_Data')
         
         % Write sum track
            write_dat(['./QA_data/'  char(model_name) '_' file_prefix '_dx_' num2str(round(dx/1E-9)) '_nm_'...
                'fx_sum_vs_t.dat'],t,fx_avg_sum_track)
            write_dat(['./QA_data/' model_name '_' file_prefix '_dx_' num2str(round(dx/1E-9)) '_nm_'...
                'fy_sum_vs_t.dat'],t,fy_avg_sum_track)
            
          % Write the profile

            write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case '_dx_' num2str(round(dx/1E-9)) '_nm_'...
                'fx_avg_vs_x.dat'],x(source1_x+2:end),fx_slab_avg(source1_x+2:end,round(Ny/2)))
            write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
                'fy_avg_vs_x.dat'],x(source1_x+2:end),fy_slab_avg(source1_x+2:end,peak_index))
            
                        write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
                'er_avg_vs_x.dat'],x(source1_x+2:end),er(source1_x+2:end,peak_index))
           
            write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
'Sx_flux_vs_t.dat'],t,Sx_flux)

write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
'Sx_flux_ref_vs_t.dat'],t,Sx_flux_ref)

          
          
       end
        
       % save entire workspace for later
       if save_mode==1
           mkdir('./saved_workspaces/')
       save_name=['./saved_workspaces/' model_name '_' file_prefix  '_sim_case_' sim_case '_dx_' num2str(round(dx/1E-9)) '_nm_Lx_' ...
           num2str(round(x_size/1E-9))...
                '_Ly_' num2str(round(y_size/dy)) '_workspace'];
            
            save(save_name)
       end
end