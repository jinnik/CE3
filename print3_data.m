function print3_data(filename, t, x, U)
    %PRINT3_DATA print 3D data to file
    %   Prints [x t U(t,x)] to a textfile named filename.
    %   Loops through x and t and prints the U value at
    %   that coordinate. Prints an empty row each time x
    %   is incremented, to make the data format suitable
    %   for pgfplots.
    fid = fopen(filename, 'w');
    for j = 1:length(x)
        for i = 1:length(t)
            fprintf(fid, '%f %f %f\r\n', x(j), t(i), U(i,j));
        end
        fprintf(fid, '\r\n');
    end
    fclose(fid);
end
