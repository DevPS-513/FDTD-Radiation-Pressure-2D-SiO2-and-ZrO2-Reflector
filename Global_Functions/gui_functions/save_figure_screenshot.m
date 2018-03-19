function [ output_args ] =save_figure_screenshot( h_call,event_data,h_axis,h_fig,filename )
% This function takes a screen shot of the axis h_axis using the
% getframe function onthe figure h_figure
%   
h_axis.Units='pixels'; % turn units into pixels, for getframe
pos=h_axis.Position; % get position in pixels of the axis to take a picture of
rect=[ 0 0 pos(3) pos(4)]; % zero the x and y index, just wude pixel width and height
pic=getframe(h_fig,pos);  % turn that frame into a pic with cdata
imwrite(pic.cdata,strcat(filename),'PNG'); % write to the filename


end

