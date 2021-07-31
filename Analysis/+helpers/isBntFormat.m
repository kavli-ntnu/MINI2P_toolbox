% Check if an input file is in BNT format
%
% Checks format of the input file. Returns TRUE if the file is in BNT
% format or FALSE otherwise. If file doesn't exist returns FALSE.
%
function result = isBntFormat(filename)
    result = false;
    if exist(filename, 'file') == 0
        return;
    end

    fid = data.safefopen(filename, 'r');
    str = fgets(fid);

    result = ischar(str) && ~isempty(strfind(lower(str), 'name')) && ~isempty(strfind(lower(str), 'version'));
end