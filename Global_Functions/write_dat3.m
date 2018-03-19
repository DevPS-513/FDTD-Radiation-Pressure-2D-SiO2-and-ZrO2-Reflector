function [  ] = write_dat3( fname,x,y,z )

delete(fname)
fileID = fopen(fname,'w+');
[nrows] = length(x);
 for row = 1:nrows
fprintf(fileID,' %d %d %d \n ' ,x(row),y(row),z(row));
 end
fclose(fileID)


end

