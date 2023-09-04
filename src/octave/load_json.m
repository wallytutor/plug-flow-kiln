function [data] = load_json(fname)
    % Load JSON data file into a structure.
    fid = fopen(fname);
    str = char(fread(fid)');
    fclose(fid);
    data = jsondecode(str);
endfunction