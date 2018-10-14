function custom_fprintf(filename, varargin)
    fid = fopen(filename, 'a');
    % strrep for nicer printing to terminal:
    fprintf(strrep(varargin{1}, ';', ' '), varargin{2:end});
    fprintf(fid, varargin{:});
    fclose(fid);
end
