function [ output_args ] = load_file_and_downsample( fname,num_points )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% load file to be downsamples
data=load(fname);

% assume an x y file, two coumns
x=data(:,1);
y=data(:,2);

N_old=length(x); % check old length

if num_points<N_old

    n_skip=round(N_old/num_points);
    n_downsample=[1:n_skip:N_old];
    
    x_new=x(n_downsample);
    y_new=y(n_downsample);
    
    write_dat(fname,x_new,y_new);
    
    
     
end


end

