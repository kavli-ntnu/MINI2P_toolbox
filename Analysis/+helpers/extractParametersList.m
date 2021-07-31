% Make list of used BNT parameters in file
%
% This function parses input file or script for entries like `<pattern>name = value`,
% and prints out the list of names. An example of such entry is: p.imageFormat = 'png'.
% It can be used to quickly check what BNT parameters are used in the script.
%
%  USAGE
%   helpers.extractParametersList(filename, pattern)
%   filename    - Path to the file to parse or script name. If script name is provided,
%                 then Matlab's function which is used to find out the corresponding
%                 file name.
%   pattern     - Begin of BNT parameters pattern. If omitted, then 'p.' is used.
%
function extractParametersList(filename, pattern)
    if nargin < 2
        pattern = 'p.';
    end
    existCode = exist(filename, 'file');
    if existCode == 0 || existCode == 2
        filename = which(filename);
    end
    parameters = containers.Map('KeyType', 'char', 'ValueType', 'char');

    fid = data.safefopen(filename, 'r');

    while ~feof(fid)
        str = fgetl(fid);
        patternPos = strfind(str, pattern);
        if ~isempty(patternPos)
            originalPos = patternPos;
            for i = 1:length(originalPos)
                patternPos = originalPos(i);

                wholeWord = regexp(str(patternPos-1:patternPos-1), '\w', 'once');
                if ~isempty(wholeWord)
                    continue;
                end
                commentsPos = strfind(str, '%');
                if ~isempty(commentsPos) && commentsPos(1) < patternPos
                    % this string is commented out, skip it
                    continue;
                end
                parameterName = regexp(str(patternPos+2:end), '\w*', 'match');
                parameterName = sprintf('%s%s', pattern, parameterName{1});

                parameters(parameterName) = parameterName;
            end
        end
    end

    cellfun(@(x) fprintf('%%  %s\n', x), values(parameters));
end