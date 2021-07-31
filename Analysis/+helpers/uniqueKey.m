% uniqueKey - Get unique key of a session.
%
% Create a MD5-based unique key based on the session basename and
% cut file(s).
%
%  USAGE
%    key = helpers.uniqueKey(basename, cuts)
%
%    basename   String. Should contain the basename.
%    cuts       Cell vector. Vector of full paths to cut file(s).
%    key        String. Resulting key. The key is independent of
%               cut names ordering, i.e. key for cuts {'1'; '2'} will be
%               the same as for cuts {'2'; '1'}.
%
% The arguments are optional. If they are not provided, then the currently
% selected session is used.
%
function key = uniqueKey(basename, cuts)
    global gCurrentTrial;
    global gBntData;

    if nargin < 2
        if gCurrentTrial <= 0 || gCurrentTrial > length(gBntData)
            error('Can not create key, because you didn''t provide input arguments and no data have been previously loaded.');
        end
        basename = gBntData{gCurrentTrial}.basename;
        cuts = gBntData{gCurrentTrial}.cuts;
    end

    [cutRows, cutCols] = size(cuts);
    cutNames = cell(1, cutRows * cutCols);
    numCuts = 0;
    for r = 1:cutRows
        for c = 1:cutCols
            if ~isempty(cuts{r, c})
                numCuts = numCuts + 1;
                [~, cutNames{numCuts}] = fileparts(cuts{r, c});
            end
        end
    end

    if numCuts > 0
        cutNames = cutNames(1:numCuts);
    end
    cutNames = sort(cutNames(:));
    opt.Format = 'HEX';
    opt.Method = 'MD5';
    opt.Input = 'ascii';
    key = DataHash(strcat(basename, cutNames{:}), opt);
end
