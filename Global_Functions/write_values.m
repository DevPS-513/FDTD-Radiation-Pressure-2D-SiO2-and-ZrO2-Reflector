function [  ] = write_values( fname,cell_data,sj,dj )
% cell data is a structure with data and its name in rows
% sj is the column number of strings

delete(fname)
fileID = fopen(fname,'w+');

[Nrows,Ncols]=size(cell_data);


j=sj;

% write cell names in top row
    for i=1:Nrows
    
    fprintf(fileID,'%s \t' ,char(cell_data(i,j)));

    end
 % write double value under each name   
fprintf(fileID,'\n');
j=dj;

    for i=1:Nrows
    
    fprintf(fileID,'%d \t' ,cell2mat(cell_data(i,j)));

    end







fclose(fileID)

end

