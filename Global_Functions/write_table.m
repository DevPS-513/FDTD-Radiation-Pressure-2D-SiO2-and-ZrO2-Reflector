function [  ] = write_table( fname,table )

delete(fname)
fileID = fopen(fname,'w+');
[nrows,ncols] = size(table);

 for row = 1:nrows     
     for col=1:ncols               
        fprintf(fileID,'%d \t' ,table(row,col));
     end 
     fprintf(fileID,' \n')
end
fclose(fileID)
end

