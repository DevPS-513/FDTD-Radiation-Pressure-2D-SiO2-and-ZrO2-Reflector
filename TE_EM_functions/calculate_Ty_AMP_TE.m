function [ Ty,t1,t2,t3,t4,t5 ] = calculate_Ty_AMP_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )


mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));



Ty=zeros(size(Ez));
Bx_avg=zeros(size(Hx));
By_avg=zeros(size(Hy));
% t1=\px(ByHx)
% we place Ty at Hx (i+1/2,j)


By_avg(i,j)=.25*(By(i,j)+By(i,j-1)+By_n_prev(i,j)+By_n_prev(i,j-1));
Bx_avg(i,j)=.25*(Bx(i,j)+Bx_n_prev(i,j)+Bx(i-1,j)+Bx_n_prev(i-1,j));


t1=-1*(1/dx)*(By_avg(i+1,j).*Bx_avg(i+1,j)-By_avg(i,j).*Bx_avg(i,j))./mu_o;


% t2=\py DzEz

t2=.5*(1/dy)*eps_o*(Ez(i,j).*Ez(i,j)-Ez(i,j-1).*Ez(i,j-1));

% t3= \py BxBx

t3=.5*(1/(2*dy))*(Bx(i,j+1).*Bx(i,j+1)-Bx(i,j-1).*Bx(i,j-1))./mu_o;

% t= \py(ByHy)

By_avg(i,j)=.5*(By(i,j)+By(i+1,j));
%Hy_avg(i,j)=.5*(Hy(i,j)+Hy(i+1,j));

t4=-0.5*(1/dy)*(By_avg(i,j).*By_avg(i,j)-By_avg(i,j-1).*By_avg(i,j-1))./mu_o;

t5=0;


%Ty(i,j)=-t1+0.5*(t2+t3-t4);

Ty(i,j)=t1+t2+t3+t4;





end

