function custom_fprintf(filename, varargin)
    %CUSTOM_FPRINTF Print to a file and terminal
    %   Prints string(s) in varargin both to a file filename
    %   and to the terminal.
    
    fid = fopen(filename, 'a');
    % strrep for nicer printing to terminal:
    fprintf(strrep(varargin{1}, ';', ' '), varargin{2:end});
    fprintf(fid, varargin{:});
    fclose(fid);
end
