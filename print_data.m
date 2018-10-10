function print_data(filename, t, U)
    fid = fopen(filename, 'w');
    for i=1:length(U)
        fprintf(fid, [repmat('%.9f ', 1, length(U(1,:))) '%.9f\n'], t(i), U(i, :));
    end
    fclose(fid);
end
