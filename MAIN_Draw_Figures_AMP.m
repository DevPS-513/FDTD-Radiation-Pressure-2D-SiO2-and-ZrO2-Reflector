close all
clear all
clc

% README: This code loads simulation data for the force on a bragg grating,
% then re-writes to text files data used for publication figures, then
% replots this data with no modification to confirm.
% MODIFICATIONS ARE:
% 1. Force values normalized to the peak AB value
% 2. z-axis normalized so that the first layer of the grating is 0
% 3. micrometers are z -axis

% Author: Max Bethune-Waddell

% 22.-23. Sept. 2016
% Tomaž Požar edited the figures for NP publication

ss=get(0,'Screensize');
sw=ss(3);

%% SI UNITS
micrometers=1E-6;
nanometers=micrometers/1000;
z_units=nanometers;

step_size='7';
sim_case='grating';
prefix='TM_prelim';

% 
% % LOAD TM fz_avg
% AB_fz_avg_TM=load(['./Figure_Data/AB_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% MN_fz_avg_TM=load(['./Figure_Data/MN_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% EL_fz_avg_TM=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% Chu_fz_avg_TM=load(['./Figure_Data/Chu_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
AMP_fz_avg_TM=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_fx_avg_vs_x.dat']);
%     EL_fz_avg_TM=load('./F_data/grating_EL_fx_avg.dat');
%     AMP_fz_avg_TM=load('./F_data/grating_AMP_fx_avg.dat');
%     Chu_fz_avg_TM=load('./F_data/grating_Chu_fx_avg.dat');
% % 
% % 
% % % LOAD TM tau peak slice
%      AB_tau_peak_TM=load(['./Figure_Data/AB_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      MN_tau_peak_TM=load(['./Figure_Data/MN_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      EL_tau_peak_TM=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      Chu_tau_peak_TM=load(['./Figure_Data/Chu_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
     AMP_tau_peak_TM=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_fy_avg_vs_x.dat']);

% load power

     TM_flux=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_Sx_flux_vs_t.dat']);

     TM_flux_ref=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_Sx_flux_ref_vs_t.dat']);

t_TM=TM_flux(:,1);
Sx_inc_TM=TM_flux(:,2);
Sx_ref_TM=TM_flux_ref(:,2);

fz_norm_TM=max(AMP_fz_avg_TM(:,2));


% % LOAD TE fz_avg
step_size='7';
sim_case='grating';
prefix='TE_prelim';

% 
% 
% AB_fz_avg_TE=load(['./Figure_Data/AB_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% MN_fz_avg_TE=load(['./Figure_Data/MN_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% EL_fz_avg_TE=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
% Chu_fz_avg_TE=load(['./Figure_Data/Chu_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fx_avg_vs_x.dat']);
AMP_fz_avg_TE=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_fx_avg_vs_x.dat']);
% 
% % LOAD TE tau peak slice
%      AB_tau_peak_TE=load(['./Figure_Data/AB_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      MN_tau_peak_TE=load(['./Figure_Data/MN_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      EL_tau_peak_TE=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
%      Chu_tau_peak_TE=load(['./Figure_Data/Chu_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
     AMP_tau_peak_TE=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_fy_avg_vs_x.dat']);

% Load TE Power

     TE_flux=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_Sx_flux_vs_t.dat']);

     TE_flux_ref=load(['./Figure_Data/AMP_' prefix '_sim_case_' sim_case '_dx_'...
step_size '_nm_Sx_flux_ref_vs_t.dat']);

t_TE=TM_flux(:,1);
Sx_inc_TE=TE_flux(:,2);
Sx_ref_TE=TE_flux_ref(:,2);

fz_norm_TE=max(AMP_fz_avg_TE(:,2));


 % DEFINE NORMALIZATION   

    dz=AMP_fz_avg_TM(2,1)-AMP_fz_avg_TM(1,1);
    source_location=1725*micrometers;
    source_location_n=round(source_location/dz);
    grating_location=2.01*micrometers;
    f_AMP_pk_TM=max(AMP_fz_avg_TM(source_location_n:end,2));
    f_AMP_pk_TE=max(AMP_fz_avg_TE(source_location_n:end,2));

% load geometry, er_1 is SiO2, er_2 is ZrO2
% units are in nm
    er_1=load('./geom_data/er_1_boxes_grating.dat');
    er_2=load('./geom_data/er_2_boxes_grating.dat');
    er_coating=load('./geom_data/er_coating_box_grating.dat');
% normalize position and magnitude
  er_2(:,1)=er_2(:,1)-dz./nanometers;
 % er_2(:,2)=er_2(:,2)./max(er_2(:,2)).*max(AB_fz_avg_TM(:,2))*2-max(AB_fz_avg_TM(:,2));
%  er_1(:,2)=er_1(:,2)./max(er_1(:,2)).*max(AB_fz_avg_TM(:,2))*2-max(AB_fz_avg_TM(:,2));
  
  
  er_2(:,2)=er_2(:,2)./max(er_2(:,2));
  er_1(:,2)=er_1(:,2)./max(er_1(:,2));
 
% Add coating to first layetr
d_coating_extra=200; % nanometers
er_1(1:2,1)=er_1(1:2,1)-d_coating_extra;
  er_1(5,1)=er_1(5,1)-d_coating_extra;

    x_limits=[ 650 4000];
% PLOT AB AND MN
f_1=figure(1);
set(f_1,'name','raw TM')

% Adjust figure position

p_1=get(f_1,'OuterPosition');
p_1(4)=2.2*p_1(4);
p_1(2)=35;
p_1(1)=0;
p_1(3)=sw/4;
set(f_1,'OuterPosition',p_1);
%% TM RESULTS
subplot(3,1,1)
hold on

h_er_1=fill(er_1(:,1),er_1(:,2),'g','Facealpha',.05);
hold on
	h_er_2=fill(er_2(:,1),er_2(:,2),'m','Facealpha',.05);
% hold on
% 	h_AB=plot(AB_fz_avg_TM(:,1)/z_units,AB_fz_avg_TM(:,2)./fz_norm_TM,'color','b');
% hold on
% 	h_MN=plot(MN_fz_avg_TM(:,1)/z_units,MN_fz_avg_TM(:,2)./fz_norm_TM,'color','r','linestyle','--');
%     
%legend([h_AB,h_MN,h_er_1,h_er_2],'AB','MN','SiO2','ZrO2','location','northeast')

%ylim([-1.5 1.5 ])
xlim(x_limits)
xlabel('z(nm)')
ylabel('\langle f_z \rangle/f^{AB}_{pk}')
title('TM')

subplot(3,1,2)

    h_er_1=fill(er_1(:,1),er_1(:,2),'g','Facealpha',.05);
hold on
	h_er_2=fill(er_2(:,1),er_2(:,2),'m','Facealpha',.05);
hold on
% 	h_EL=plot(EL_fz_avg_TM(:,1)/z_units,EL_fz_avg_TM(:,2)./fz_norm_TM,'black','linewidth',1.5);
% hold on
% 	h_Chu=plot(Chu_fz_avg_TM(:,1)/z_units,Chu_fz_avg_TM(:,2)./fz_norm_TM,'color','magenta','linestyle','-.');
% hold on
	h_AMP=plot(AMP_fz_avg_TM(:,1)/z_units,AMP_fz_avg_TM(:,2)./fz_norm_TM,'color','green','linestyle','--');
    
%legend([h_EL,h_Chu,h_AMP],'EL','Chu','AMP','location','northeast')
% ylim([-1.5 1.5 ])
xlim(x_limits)
xlabel('z(nm)')
ylabel('\langle f_z \rangle/f^{AB}_{pk}')

subplot(3,1,3)
er_fac=1./max(AMP_fz_avg_TM(:,2))*max(AMP_tau_peak_TM(:,2));
    h_er_1=fill(er_1(:,1),er_1(:,2)*er_fac,'g','Facealpha',.05);
    hold on
	h_er_2=fill(er_2(:,1),er_2(:,2)*er_fac,'m','Facealpha',.05);
hold on
%  	h_AB=plot(AB_tau_peak_TM(:,1)/z_units,AB_tau_peak_TM(:,2)./fz_norm_TM,'b');
%   hold on
%  	h_AB=plot(MN_tau_peak_TM(:,1)/z_units,AB_tau_peak_TM(:,2)./fz_norm_TM,'--r');
% hold on
% 	h_EL=plot(EL_tau_peak_TM(:,1)/z_units,EL_tau_peak_TM(:,2)./fz_norm_TM,'black','linewidth',1.5);
% hold on
	h_AMP=plot(AMP_tau_peak_TM(:,1)/z_units,AMP_tau_peak_TM(:,2)./fz_norm_TM,'g');
% hold on
% 	h_Chu=plot(Chu_tau_peak_TM(:,1)/z_units,Chu_tau_peak_TM(:,2)./fz_norm_TM,'--magenta');
%legend([h_AB,h_MN,h_EL,h_Chu,h_AMP],'AB & MN','EL','Chu','AMP','location','northeast')

%     ylim([-1.5 1.5 ])
xlim(x_limits)
    xlabel('z(nm)')
    ylabel('\langle f_r \rangle/f^{AB}_{pk}')


f_2=figure(2);
set(f_2,'name','raw TE')
p_2=get(f_2,'OuterPosition');
% make second figure to the right of first figure
p_2(1)=p_1(1)+p_1(3);
p_2(2)=p_1(2);
p_2(3)=p_1(3);
p_2(4)=p_1(4);
set(f_2,'OuterPosition',p_2);
%% TE RESULTS
subplot(3,1,1)
hold on

h_er_1=fill(er_1(:,1),er_1(:,2),'g','Facealpha',.05);
hold on
	h_er_2=fill(er_2(:,1),er_2(:,2),'m','Facealpha',.05);
% hold on
% 	h_AB=plot(AB_fz_avg_TE(:,1)/z_units,AB_fz_avg_TE(:,2)./fz_norm_TE,'color','b');
% hold on
% 	h_MN=plot(MN_fz_avg_TE(:,1)/z_units,MN_fz_avg_TE(:,2)./fz_norm_TE,'color','r','linestyle','--');
%     
%legend([h_AB,h_MN,h_er_1,h_er_2],'AB','MN','SiO2','ZrO2','location','northeast')

%ylim([-1.5 1.5 ])
xlim(x_limits)
xlabel('z(nm)')
ylabel('\langle f_z \rangle/f^{AB}_{pk}')
title('TE')

subplot(3,1,2)

    h_er_1=fill(er_1(:,1),er_1(:,2),'g','Facealpha',.05);
hold on
	h_er_2=fill(er_2(:,1),er_2(:,2),'m','Facealpha',.05);
% hold on
% 	h_EL=plot(EL_fz_avg_TE(:,1)/z_units,EL_fz_avg_TE(:,2)./fz_norm_TE,'black','linewidth',1.5);
% hold on
% 	h_Chu=plot(Chu_fz_avg_TE(:,1)/z_units,Chu_fz_avg_TE(:,2)./fz_norm_TE,'color','magenta','linestyle','-.');
% hold on
	h_AMP=plot(AMP_fz_avg_TE(:,1)/z_units,AMP_fz_avg_TE(:,2)./fz_norm_TE,'color','green','linestyle','--');
    
%legend([h_EL,h_Chu,h_AMP],'EL','Chu','AMP','location','northeast')
% ylim([-1.5 1.5 ])
xlim(x_limits)
xlabel('z(nm)')
ylabel('\langle f_z \rangle/f^{AB}_{pk}')

subplot(3,1,3)
%er_fac=1./max(AB_fz_avg_TM(:,2))*max(EL_tau_peak_TE(:,2));
  %  h_er_1=fill(er_1(:,1),er_1(:,2)*er_fac,'g','Facealpha',.05);
  %  hold on
  %  h_er_2=fill(er_2(:,1),er_2(:,2)*er_fac,'m','Facealpha',.05);
% hold on
% 	h_AB=plot(AB_tau_peak_TE(:,1)/z_units,AB_tau_peak_TE(:,2)./fz_norm_TE,'b');
% % hold on
%  	h_MN=plot(MN_tau_peak_TE(:,1)/z_units,MN_tau_peak_TE(:,2)./fz_norm_TE,'--r');
% hold on
% 	h_EL=plot(EL_tau_peak_TE(:,1)/z_units,EL_tau_peak_TE(:,2)./fz_norm_TE,'.-black');
% hold on
	h_AMP=plot(AMP_tau_peak_TE(:,1)/z_units,AMP_tau_peak_TE(:,2)./fz_norm_TE,'g');
% hold on
% 	h_Chu=plot(Chu_tau_peak_TE(:,1)/z_units,Chu_tau_peak_TE(:,2)./fz_norm_TE,'--magenta');
%legend([h_AB,h_MN,h_EL,h_Chu,h_AMP],'AB','MN','EL','Chu','AMP','location','northeast')

%     ylim([-1.5 1.5 ])
xlim(x_limits)
    xlabel('z(nm)')
    ylabel('\langle f_r \rangle/f^{AB}_{pk}')
    
    % plot power
%   
%     t_TE=TM_flux(:,1);
% Sx_inc_TE=TE_flux(:,2);
% Sx_ref_TE=TE_flux_ref(:,2);
    
figure(3)
plot(t_TE,Sx_inc_TE,'b')
% hold on
% plot(t_TE,Sx_ref_TE,'r')
hold on
plot(t_TE,Sx_inc_TM,'--black')
% hold on
% plot(t_TE,Sx_ref_TM,'.g')
legend('TE','TM')
   

% Plot total force 

% Write sum track
%            write_dat(['./QA_data/'  char(model_name) '_' file_prefix '_dx_' num2str(round(dx/1E-9)) '_nm_'...
%                 'fx_sum_vs_t.dat'],t,fx_avg_sum_track)
 %            write_dat(['./QA_data/' model_name '_' file_prefix '_dx_' num2str(round(dx/1E-9)) '_nm_'...
  %               'fy_sum_vs_t.dat'],t,fy_avg_sum_track)
   %          
          % Write the profile

     %        write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case '_dx_' num2str(round(dx/1E-9)) '_nm_'...
       %          'fx_avg_vs_x.dat'],x(source1_x+2:end),fx_slab_avg(source1_x+2:end,round(Ny/2)))
         %    write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
           %      'fy_avg_vs_x.dat'],x(source1_x+2:end),fy_slab_avg(source1_x+2:end,peak_index))
            
             %            write_dat(['./Figure_Data/' model_name '_' file_prefix '_sim_case_' sim_case  '_dx_' num2str(round(dx/1E-9)) '_nm_'...
         %        'er_avg_vs_x.dat'],x(source1_x+2:end),er(source1_x+2:end,peak_index))
           
%load TM data	


% 	step_size='8';
% 	sim_case='grating';
% 	prefix='TM';
% 
% 
% 	fx_sum_TM=load(	['./QA_Data/AB_' prefix '_dx_'...
% 	step_size '_nm_fx_sum_vs_t.dat'])   ;
% 
% 	t_TM=fx_sum_TM(:,1);
% 	fx_sum_TM=fx_sum_TM(:,2);
% 
% 
% 	fy_sum_TM=load(	['./QA_Data/AB_' prefix '_dx_'...
% 	step_size '_nm_fy_sum_vs_t.dat'])   ;
% 
% 
% 	%t_TM=fx_sum_TM(:,1);
% 	fy_sum_TM=fy_sum_TM(:,2);
% 	
% 	% load TE data
% 
% 	
% 	step_size='8';
% 	sim_case='grating';
% 	prefix='TE';
% 
% 
% 	fx_sum_TE=load(	['./QA_Data/AB_' prefix '_dx_'...
% 	step_size '_nm_fx_sum_vs_t.dat'])   ;
% 
% 
% 	t_TE=fx_sum_TE(:,1);
% 	fx_sum_TE=fx_sum_TE(:,2);
% 
% 
% 	fy_sum_TE=load(	['./QA_Data/AB_' prefix '_dx_'...
% 	step_size '_nm_fy_sum_vs_t.dat'])   ;
% 
% 
% 	%t_TM=fx_sum_TM(:,1);
% 	fy_sum_TE=fy_sum_TE(:,2);
% 		   
% 		   
% 		   
% 
%     
%     
%     
% fig_4=figure(4);
% subplot(2,1,1)
% 	plot(t_TM,fx_sum_TM,'blue')
% 	hold on 
% 	plot(t_TE,fx_sum_TE,'.r')
% 	ylabel('fx')
% 	xlabel('t')
% 
% legend('11 TM','11 TE')
% subplot(2,1,2)
% 	plot(t_TM,fy_sum_TM,'blue')
% 	hold on
% 	plot(t_TE,fy_sum_TE,'.r')
% 	legend('11 TM','11 TE')
% 	ylabel('fy')
% 	xlabel('t')
%  
%     
% %% Load distributions at different step sizes
% 
% % LOAD TE as a reference
%     step_size='12';
%     sim_case='grating';
%     prefix='TE';
% 
%     AB_fy_12_TE=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
%     step_size '_nm_fy_avg_vs_x.dat']);
% 
%     z_AB_fy_12_TE=AB_fy_12_TE(:,1);
%     AB_fy_12_TE=AB_fy_12_TE(:,2);   
% 
% 
%     
% step_size='8';
% sim_case='grating';
% prefix='TM';
% 
%   
%     
% AB_fy_8=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']);
% 
% z_AB_fy_8=AB_fy_8(:,1);
% 
% AB_fy_8=AB_fy_8(:,2);
% 
% 
% step_size='12';
% sim_case='grating';
% prefix='TM';
%     
% AB_fy_12=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']); 
% 
% z_AB_fy_12=AB_fy_12(:,1);
% AB_fy_12=AB_fy_12(:,2);
% 
% 
% step_size='14';
% sim_case='grating';
% prefix='TM';
%     
% AB_fy_14=load(['./Figure_Data/EL_' prefix '_sim_case_' sim_case '_dx_'...
% step_size '_nm_fy_avg_vs_x.dat']); 
%     
%     
% z_AB_fy_14=AB_fy_14(:,1);
% AB_fy_14=AB_fy_14(:,2);
% 
% fig_5=figure(5);
% h_8=plot(z_AB_fy_8,AB_fy_8,'r');
% % hold on
% % h_12=plot(z_AB_fy_12,AB_fy_12,'blue');
% hold on
% h_14=plot(z_AB_fy_14,AB_fy_14,'blue');
% hold on
% h_12_TE=plot(z_AB_fy_12_TE,AB_fy_12_TE,'--black');
% 
% legend([h_8,h_14,h_12_TE],'8 nm TM','14 nm TM','12 nm TE')
% 
% xlabel('z (m)')
% ylabel('f_z [N] (a.u sim)')
% xlim([.6E-6 2.5E-6])
% %     
% % % MODIFY FILES TO ERASE FORCE FROM SOURCE AND PML
% % erase_source_force_TM=ones(size(AB_fz_avg_TM(:,2)));
% % erase_source_force_TM(1:source_location_n)=0;
% % erase_source_force_TE=ones(size(AB_fz_avg_TE(:,2)));
% % erase_source_force_TE(1:source_location_n)=0;
% % % mkdir('./Paper_Data');
% % 
% % 
% % % WRITE fz of GRATING TM
% % write_dat('./Paper_Data/grating_AB_fz_avg_vs_mum_TM.dat',...
% %         AB_fz_avg_TM(:,1)/z_units-grating_location,...
% %         erase_source_force_TM.*AB_fz_avg_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_MN_fz_avg_vs_mum_TM.dat',...
% %         MN_fz_avg_TM(:,1)/z_units-grating_location,...
% %         erase_source_force_TM.*MN_fz_avg_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_EL_fz_avg_vs_mum_TM.dat',...
% %         EL_fz_avg_TM(:,1)/z_units-grating_location,...
% %         erase_source_force_TM.*EL_fz_avg_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_Chu_fz_avg_vs_mum_TM.dat',...
% %         Chu_fz_avg_TM(:,1)/z_units-grating_location,...
% %         erase_source_force_TM.*Chu_fz_avg_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_AMP_fz_avg_vs_mum_TM.dat',...
% %         AMP_fz_avg_TM(:,1)/z_units-grating_location,...        
% %         erase_source_force_TM.*AMP_fz_avg_TM(:,2)./f_AB_pk_TM);
% %     
% % % WRITE fz of GRATING TE
% % write_dat('./Paper_Data/grating_AB_fz_avg_vs_mum_TE.dat',...
% %         AB_fz_avg_TE(:,1)/z_units-grating_location,...
% %         erase_source_force_TE.*AB_fz_avg_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_MN_fz_avg_vs_mum_TE.dat',...
% %         MN_fz_avg_TE(:,1)/z_units-grating_location,...
% %         erase_source_force_TE.*MN_fz_avg_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_EL_fz_avg_vs_mum_TE.dat',...
% %         EL_fz_avg_TE(:,1)/z_units-grating_location,...
% %         erase_source_force_TE.*EL_fz_avg_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_Chu_fz_avg_vs_mum_TE.dat',...
% %         Chu_fz_avg_TE(:,1)/z_units-grating_location,...
% %         erase_source_force_TE.*Chu_fz_avg_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_AMP_fz_avg_vs_mum_TE.dat',...
% %         AMP_fz_avg_TE(:,1)/z_units-grating_location,...        
% %         erase_source_force_TE.*AMP_fz_avg_TE(:,2)./f_AB_pk_TE);	
% % % WRITE PEAK TAU SLICE TM
% %     write_dat('./Paper_Data/grating_AB_tau_pk_vs_mum_TM.dat',...
% %         AB_tau_peak_TM(:,1)/z_units-grating_location,AB_tau_peak_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_MN_tau_pk_vs_mum_TM.dat',...
% %         MN_tau_peak_TM(:,1)/z_units-grating_location,MN_tau_peak_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_EL_tau_pk_vs_mum_TM.dat',...
% %         EL_tau_peak_TM(:,1)/z_units-grating_location,EL_tau_peak_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_Chu_tau_pk_vs_mum_TM.dat',...
% %         Chu_tau_peak_TM(:,1)/z_units-grating_location,Chu_tau_peak_TM(:,2)./f_AB_pk_TM);
% %     write_dat('./Paper_Data/grating_AMP_tau_pk_vs_mum_TM.dat',...
% %         Chu_tau_peak_TM(:,1)/z_units-grating_location,Chu_tau_peak_TM(:,2)./f_AB_pk_TM);
% % % WRITE PEAK TAU SLICE TE
% % write_dat('./Paper_Data/grating_AB_tau_pk_vs_mum_TE.dat',...
% %         AB_tau_peak_TE(:,1)/z_units-grating_location,AB_tau_peak_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_MN_tau_pk_vs_mum_TE.dat',...
% %         MN_tau_peak_TE(:,1)/z_units-grating_location,MN_tau_peak_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_EL_tau_pk_vs_mum_TE.dat',...
% %         EL_tau_peak_TE(:,1)/z_units-grating_location,EL_tau_peak_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_Chu_tau_pk_vs_mum_TE.dat',...
% %         Chu_tau_peak_TE(:,1)/z_units-grating_location,Chu_tau_peak_TE(:,2)./f_AB_pk_TE);
% %     write_dat('./Paper_Data/grating_AMP_tau_pk_vs_mum_TE.dat',...
% %         AMP_tau_peak_TE(:,1)/z_units-grating_location,AMP_tau_peak_TE(:,2)./f_AB_pk_TE);
% % 
% % 		
% % 		
% % 		
% % % Re-write geometry with better names
% % write_dat('./Paper_Data/er_SiO2.dat',...
% %         er_1(:,1)/z_units-grating_location,er_1(:,2));5
% % write_dat('./Paper_Data/er_ZrO2.dat',...
% %         er_2(:,1)/z_units-grating_location,er_2(:,2));
% % write_dat('./Paper_Data/er_SiO2_coating.dat',...
% %         er_coating(:,1)/z_units-grating_location,er_coating(:,2));		
% % 
% % %% RE-LOAD DATA AND PLOT
% % % Load fz, note this is "fx" in the data file name
% % AB_fz_avg_TM=load('./Paper_Data/grating_AB_fz_avg_vs_mum_TM.dat');
% % MN_fz_avg_TM=load('./Paper_Data/grating_MN_fz_avg_vs_mum_TM.dat');
% % EL_fz_avg_TM=load('./Paper_Data/grating_EL_fz_avg_vs_mum_TM.dat');
% % AMP_fz_avg_TM=load('./Paper_Data/grating_AMP_fz_avg_vs_mum_TM.dat');
% % Chu_fz_avg_TM=load('./Paper_Data/grating_Chu_fz_avg_vs_mum_TM.dat');
% % 
% % % load tau peak slice
% % AB_tau_peak_TM=load('./Paper_Data/grating_AB_tau_pk_vs_mum_TM.dat');
% % MN_tau_peak_TM=load('./Paper_Data/grating_MN_tau_pk_vs_mum_TM.dat');
% % EL_tau_peak_TM=load('./Paper_Data/grating_EL_tau_pk_vs_mum_TM.dat');
% % Chu_tau_peak_TM=load('./Paper_Data/grating_Chu_tau_pk_vs_mum_TM.dat');
% % AMP_tau_peak_TM=load('./Paper_Data/grating_AMP_tau_pk_vs_mum_TM.dat');
% % 
% % %% RE-LOAD DATA AND PLOT
% % % Load fz, note this is "fx" in the data file name
% % AB_fz_avg_TE=load('./Paper_Data/grating_AB_fz_avg_vs_mum_TE.dat');
% % MN_fz_avg_TE=load('./Paper_Data/grating_MN_fz_avg_vs_mum_TE.dat');
% % EL_fz_avg_TE=load('./Paper_Data/grating_EL_fz_avg_vs_mum_TE.dat');
% % AMP_fz_avg_TE=load('./Paper_Data/grating_AMP_fz_avg_vs_mum_TE.dat');
% Chu_fz_avg_TE=load('./Paper_Data/grating_Chu_fz_avg_vs_mum_TE.dat');
% 
% % load tau peak slice
% AB_tau_peak_TE=load('./Paper_Data/grating_AB_tau_pk_vs_mum_TE.dat');
% MN_tau_peak_TE=load('./Paper_Data/grating_MN_tau_pk_vs_mum_TE.dat');
% EL_tau_peak_TE=load('./Paper_Data/grating_EL_tau_pk_vs_mum_TE.dat');
% Chu_tau_peak_TE=load('./Paper_Data/grating_Chu_tau_pk_vs_mum_TE.dat');
% AMP_tau_peak_TE=load('./Paper_Data/grating_AMP_tau_pk_vs_mum_TE.dat');
% 
% 
% % Load Geometry
% er_SiO2=load('./Paper_data/er_SiO2.dat');
% er_SiO2_coating=load('./Paper_data/er_SiO2_coating.dat');
% er_ZrO2=load('./Paper_data/er_ZrO2.dat');
% 
% % -------------------------------------------------------------------------
% % NP Plot (Edited by Tomaž Požar)
% % -------------------------------------------------------------------------
% 
% 
% % Find peaks
% 
% [max_vals peak_locs]=findpeaks(AB_fz_avg_TM(:,2));
% [min_vals min_locs]=findpeaks(-1*AB_fz_avg_TM(:,2));
% 
% peak_locs=peak_locs(max_vals>.003);
% min_locs=min_locs(min_vals<-.002);
% 
% [global_max, global_loc]=max(AB_fz_avg_TM(:,2));
% 
% 
% 
% peak_diff=diff(AB_fz_avg_TM(peak_locs,1));
% 
% layer_width=.345; % found from peaks
% 
% er_new=[];
% % 
% x_er_new_SiO2=[er_SiO2_coating(:,1)'];
%  y_er_new_SiO2=[er_SiO2_coating(:,2)'];
%  
%  x_er_new_ZrO2=[];
%  y_er_new_ZrO2=[];
% 
%  peak_locs=[global_loc peak_locs'];
%  
% for j=5:length(peak_locs)
%     
%     x1=AB_fz_avg_TM(peak_locs(j),1);
%     
%     x2=x1+layer_width;
%     
%     x_box_SiO2=[x1 x1 x2 x2 x1];
%     y_box_SiO2=[-1 1 1 -1 -1];
%     
%     % now go backwards for other coating 136 nm
%     x_box_other=[x1 x1 x2-.136 x2-.136 x1];
%     
% %     if j>6
% %     x_er_new_ZrO2=[x_er_new_ZrO2 x_box_other];
% %     y_er_new_ZrO2=[y_er_new_ZrO2 y_box_SiO2];
% %     end
%    
%     
% end
% 
%  x_er_new_SiO2(1,3:4)=0.885;
%  
%  
% 
% figure(77)
% plot(AB_fz_avg_TM(:,1),AB_fz_avg_TM(:,2))
% hold on
% fill(x_er_new_SiO2,y_er_new_SiO2,'g','facealpha',.1)
% hold on
% fill(x_er_new_SiO2,y_er_new_SiO2,'m','facealpha',.1)
% %interface_locs=mean
% 
% 
% 
% f_3=figure(3);
% p_3=get(f_3,'OuterPosition');
% % make second figure to the right of first figure
% p_3(1)=p_2(1)+p_2(3);
% p_3(2)=p_2(2);
% p_3(3)=p_2(3);
% p_3(4)=p_2(4);
% set(f_3,'OuterPosition',p_3);
% 
% force_x_offset=(.0553-.054)*0;
% 
% subplot(2,1,1)
% fill(er_ZrO2(:,1)+.011,er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','none','Edgealpha',1);
% hold on
% fill(er_ZrO2(:,1),er_ZrO2(:,2),'white','Facealpha',1,'EdgeColor','none','Edgealpha',1);
% hold on
% 
% 	h_er_1=fill(er_SiO2(:,1),er_SiO2(:,2),'g','Facealpha',.1,'EdgeColor','None');
% hold on
% 	h_er_2=fill(er_ZrO2(:,1),er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','None');
% hold on
%     h_er_c=fill(er_SiO2_coating(:,1),er_SiO2_coating(:,2),'g','Facealpha',.1,'EdgeColor','None');
% 
% hold on
% 	h_AB = plot(AB_fz_avg_TM(:,1)+force_x_offset,AB_fz_avg_TM(:,2),'k'); % AB and MN
% hold on
% 	h_EL = plot(EL_fz_avg_TM(:,1)+force_x_offset,EL_fz_avg_TM(:,2),'b');  % EL
%     hold on
% 	h_Chu = plot(AMP_fz_avg_TM(:,1)+force_x_offset,Chu_fz_avg_TM(:,2),'--m');  % Chu and AMP
% hold on
% 	h_AMP = plot(AMP_fz_avg_TM(:,1)+force_x_offset,AMP_fz_avg_TM(:,2),'--r');  % Chu and AMP
% 
% legend([h_EL,h_AB,h_AMP],'EL','AB and MN','Chu and AMP','location','northeast')
% ylim([-0.5 1])
% xlim([-0.5 4.2])
% set(gca,'Xtick',0:4);
% xlabel('z(\mu{m})')
% ylabel('\langle f_z \rangle/f^{AB}_{z,pk}')
% title('TM')
% 
% 
% subplot(2,1,2)
% fill(er_ZrO2(:,1)+.011,er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','none','Edgealpha',1);
% hold on
% fill(er_ZrO2(:,1),er_ZrO2(:,2),'white','Facealpha',1,'EdgeColor','none','Edgealpha',1);
% hold on
% 
%     h_er_1=fill(er_SiO2(:,1),er_SiO2(:,2),'g','Facealpha',.1,'EdgeColor','None');
% hold on
% 	h_er_2=fill(er_ZrO2(:,1),er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','None');
% hold on
%     h_er_c=fill(er_SiO2_coating(:,1),er_SiO2_coating(:,2),'g','Facealpha',.1,'EdgeColor','None');
% hold on    
% 	h_AB = plot(AB_tau_peak_TM(:,1)+force_x_offset,AB_tau_peak_TM(:,2),'k'); % AB and MN
% hold on
% 	h_EL = plot(EL_tau_peak_TM(:,1)+force_x_offset,EL_tau_peak_TM(:,2),'b'); % EL
% hold on
%     h_Chu = plot(Chu_tau_peak_TM(:,1)+force_x_offset,Chu_tau_peak_TM(:,2),'--m'); % Chu and AMP
%     hold on
% 	h_AMP = plot(AMP_tau_peak_TM(:,1)+force_x_offset,AMP_tau_peak_TM(:,2),'--r'); % Chu and AMP
% 
% 
% legend([h_EL,h_AB,h_AMP],'EL','AB and MN','Chu and AMP','location','northeast')
% ylim([-0.5 1])
% xlim([-0.5 4.2])
% set(gca,'Xtick',0:4)
% xlabel('z(\mu{m})')
% ylabel('\langle f_r \rangle/f^{AB}_{z,pk}')
% title('TM')
% AB_total=sum(AB_tau_peak_TM(:,2)./max(abs(AB_tau_peak_TM(:,2))));
% MN_total=sum(MN_tau_peak_TM(:,2)./max(abs(MN_tau_peak_TM(:,2))));
% EL_total=sum(EL_tau_peak_TM(:,2)./max(abs(EL_tau_peak_TM(:,2))));
% AMP_total=sum(AMP_tau_peak_TM(:,2)./max(abs(AMP_tau_peak_TM(:,2))));
% 
% disp('relative to max(tau)')
% disp([ ' sum of AB tau peak :' num2str(AB_total) ]);
% disp([ ' sum of MN tau peak :' num2str(MN_total) ]);
% disp([ ' sum of EL tau peak :' num2str(EL_total) ]);
% disp([ ' sum of AMP tau peak :' num2str(AMP_total) ]);
% 
% 
% disp('sum(tau) relative to sum(fz)')
% AB_total=sum(AB_tau_peak_TM(:,2))./sum((AB_fz_avg_TM(:,2)));
% MN_total=sum(MN_tau_peak_TM(:,2))./sum((MN_fz_avg_TM(:,2)));
% EL_total=sum(EL_tau_peak_TM(:,2))./sum((EL_fz_avg_TM(:,2)));
% AMP_total=sum(AMP_tau_peak_TM(:,2))./sum((AMP_fz_avg_TM(:,2)));
% 
% disp([ ' sum of AB tau peak :' num2str(AB_total) ]);
% disp([ ' sum of MN tau peak :' num2str(MN_total) ]);
% disp([ ' sum of EL tau peak :' num2str(EL_total) ]);
% disp([ ' sum of AMP tau peak :' num2str(AMP_total) ]);
% 
% 
% 
% 
% f_4=figure(4);
% p_4=get(f_4,'OuterPosition');
% % make second figure to the right of first figure
% p_4(1)=p_3(1)+p_3(3);
% p_4(2)=p_3(2);
% p_4(3)=p_3(3);
% p_4(4)=p_3(4);
% set(f_4,'OuterPosition',p_4);
% 
% 
% subplot(2,1,1)
% % Cover up plotting offset
% fill(er_ZrO2(:,1)+.011,er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','none','Edgealpha',1);
% hold on
% fill(er_ZrO2(:,1),er_ZrO2(:,2),'white','Facealpha',1,'EdgeColor','none','Edgealpha',1);
% hold on
% 
% 	h_er_1=fill(er_SiO2(:,1),er_SiO2(:,2),'g','Facealpha',.1,'EdgeColor','None');
% hold on
% 	h_er_2=fill(er_ZrO2(:,1),er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','None');
% hold on
%     h_er_c=fill(er_SiO2_coating(:,1),er_SiO2_coating(:,2),'g','Facealpha',.1,'EdgeColor','None');
% 
% hold on
% 	h_AB = plot(AB_fz_avg_TE(:,1),AB_fz_avg_TE(:,2),'k'); % AB and MN
% hold on
% 	h_EL = plot(EL_fz_avg_TE(:,1),EL_fz_avg_TE(:,2),'b');  % EL
% hold on
% 	h_Chu = plot(Chu_fz_avg_TE(:,1),Chu_fz_avg_TE(:,2),'--m');  % Chu and AMP
% hold on
% 	h_AMP = plot(AMP_fz_avg_TE(:,1),AMP_fz_avg_TE(:,2),'--r');  % Chu and AMP
% 
% 
% legend([h_EL,h_AB,h_AMP],'EL','AB and MN','Chu and AMP','location','northeast')
% ylim([-0.5 1])
% xlim([-0.5 4.2])
% set(gca,'Xtick',0:4);
% xlabel('z(\mu{m})')
% ylabel('\langle f_z \rangle/f^{AB}_{z,pk}')
% title('TE')
% 
% subplot(2,1,2)
% 
% fill(er_ZrO2(:,1)+.011,er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','none','Edgealpha',1);
% hold on
% fill(er_ZrO2(:,1),er_ZrO2(:,2),'white','Facealpha',1,'EdgeColor','none','Edgealpha',1);
% hold on
%     h_er_1=fill(er_SiO2(:,1),er_SiO2(:,2),'g','Facealpha',.1,'EdgeColor','None');
% hold on
% 	h_er_2=fill(er_ZrO2(:,1),er_ZrO2(:,2),'m','Facealpha',.1,'EdgeColor','None');
% hold on
%     h_er_c=fill(er_SiO2_coating(:,1),er_SiO2_coating(:,2),'g','Facealpha',.1,'EdgeColor','None');
% 
% hold on    
% 	h_AB = plot(AB_tau_peak_TE(:,1),AB_tau_peak_TE(:,2),'k'); % AB and MN
% hold on
% 	h_EL = plot(EL_tau_peak_TE(:,1),EL_tau_peak_TE(:,2),'black.'); % EL
% hold on
% 	h_AMP = plot(AMP_tau_peak_TE(:,1),AMP_tau_peak_TE(:,2),'--r'); % Chu and AMP
% %hold on
% 	%h_Chu=plot(Chu_tau_peak_TE(:,1),Chu_tau_peak_TE(:,2),'magenta');    
% 
% legend([h_EL,h_AB,h_AMP,h_Chu],'EL','AB and MN','Chu and AMP','location','northeast')
% ylim([-0.5 1])
% xlim([-0.5 4.2])
% set(gca,'Xtick',0:4)
% xlabel('z(\mu{m})')
% ylabel('\langle f_r \rangle/f^{AB}_{z,pk}')
% title('TE')
% 
% 
% figure
% 
% EL_Fy_Data=load('./F_data/AB_Fy_slab_avg_grating_TE_cw.dat');
% 
% x_Fy=EL_Fy_Data(:,1);
% Fy_EL=EL_Fy_Data(:,2);
% 
% plot(x_Fy,Fy_EL)
% 
% 
