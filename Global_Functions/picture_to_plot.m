function [ x,y ] = picture_to_plot(name,ext,xmin,xmax,ymin,ymax,RGB_mat,N_points )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% FILENAME

 X=imread(strcat(name,'.',ext));
 
 % set which RGB values to find
% R_find=51;
% G_find=102;
% B_find=204;

R_find=RGB_mat(1);
G_find=RGB_mat(2);
B_find=RGB_mat(3);

% x and y limits of plot
x_1=xmin;     % deg
x_end=xmax;   % deg

y_1=ymin;       % pN
y_end=ymax;     % pN


 
%  % Display image in matlab
%  f_1=figure(1);
%  subplot(1,3,1)
%  imshow(X);
%  axis xy
%  title(' original image, matlab flips it upside down...')
%  set(f_1,'OuterPosition',[0 35 1000 600])


 % EXTRACT PIXELS
R=X(:,:,1);
G=X(:,:,2);
B=X(:,:,3);
W=1.*(1.*(R==R_find).*1.*(G==G_find).*1.*(B==B_find));
%% COMMENT THESE THREE LINES IF YOUR PLOT IS UPSIDE DOWN ON THE OUTPUT!
W=fliplr(W);
W=flipud(W);
W=fliplr(W);

% check if your image needs to be flipped around..(this one does, remedy
% above)
% plot ones where the desired pixels are
% f_1=figure(1);
% subplot(1,3,2)
% surf(W);
% shading flat
% view([ 0 90])
% title('image data')


% define axis
[Ny,Nx]=size(W);
dx=(x_end-x_1)./(Nx-1);
x=[0:1:Nx-1]*dx+x_1;
dy=(y_end-y_1)./(Ny-1);
y=[0:1:Ny-1]*dy+y_1;


% put locations into a matrix
[J,I,Q]=find(W);

% fname='results.txt';
% delete(fname)
% fileID = fopen(fname,'w+');
% [nrows] = length((I));
%  for row = 1:nrows
% fprintf(fileID,' %d %d \n ' ,x(I(row)),y(J(row)));
%  end
% fclose(fileID)

% plot data

% 
% figure(1)
% subplot(1,3,3)
% plot(x(I),y(J),'color','r')
% xlabel('sep (in)')
% ylabel('force (lb)')
% title('interpolated plot!')

x=x(I);
y=y(J);

[x_new IA ,IC]=unique(x);
y_new=y(IA);

dx_new=(x_new(end)-x_new(1))/(N_points);
xq=x_new(1):dx_new:x_new(end);

yq=interp1(x_new,y_new,xq,'linear','extrap');

x=xq;
y=yq;



end

