% Modified Matlab function fileparts that is able to detect dots in file names
%
% For the input 'c:\path\complex.file and long name' Matlab's fileparts function
% will return 'complex' as name and '.file and long name' as ext. This modified
% version will in general return 'complex.file and long name' as name and empty
% string as ext.
%
% For see also fileparts
%   
function [pathstr, name, ext] = fileparts(filename)
    [pathstr, name, ext] = fileparts(filename);
    if length(ext) > 4
        name = [name ext];
        ext = '';
    end
end