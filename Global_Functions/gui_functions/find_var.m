function [ val,loc ] = find_var( struct_mat,data,name,name_index )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

loc=find(strcmp(struct_mat(:,name_index),name));

val= data(loc); 

if(iscell(val))
    
   val=cell2mat(val); 
end
end

