function print_data(filename, t, U)
    %PRINT_DATA print 2D data to file
    %   Prints [t U] to a textfile named filename.
    %   U can have 1 or more columns.
    fid = fopen(filename, 'w');
    for i=1:length(U)
        fprintf(fid, [repmat('%.9f ', 1, length(U(1,:))) '%.9f\n'], t(i), U(i, :));
    end
    fclose(fid);
end
