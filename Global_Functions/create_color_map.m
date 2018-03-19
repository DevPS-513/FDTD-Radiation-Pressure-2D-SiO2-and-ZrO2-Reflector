function [ map ] = create_color_map( negative_color,middle_color,positive_color,n_colors )
% input values assumed to be [R G B], out of 1
% negative color,   the color at the maximum negative value
% middle color,     the color at the zero point
% positive color,   the color at the maximum positive value


% specify "reference matricies"
% to turn three colors into an n_color colormap
R_ref=[negative_color(1) middle_color(1) positive_color(1)];
G_ref=[negative_color(2) middle_color(2) positive_color(2)];
B_ref=[negative_color(3) middle_color(3) positive_color(3)];


% specify x axis for reference amtricies   
x_ref=[1 round(n_colors/2) n_colors];
    
R=interp1(x_ref,R_ref,[1:n_colors]);
G=interp1(x_ref,G_ref,[1:n_colors]);
B=interp1(x_ref,B_ref,[1:n_colors]);

map=[R;G;B;]';

end

