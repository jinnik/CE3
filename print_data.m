function print_data(filename, t, x, U)
fid = fopen(filename, 'w');
for j = 1:length(x)
    for i = 1:length(t)
        fprintf(fid, '%f %f %f\r\n', x(j), t(i), U(i,j)); %, i, j, U(i,j));
    end
    fprintf(fid, '\r\n');
end
fclose(fid);
end
