function [data] = dump_json(fname, data)
    % Dump structure into a JSON data file.
    fid = fopen(fname, "w");
    fputs(fid, jsonencode(data));
    fclose(fid);
endfunction