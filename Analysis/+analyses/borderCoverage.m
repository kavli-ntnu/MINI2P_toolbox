% Calculate border coverage for detected fields
%
% This function calculates firing map border coverage that is further
% used in calculation of a border score.
%
%  USAGE
%   coverage = analyses.borderCoverage(fields, <options>)
%   fields          Array of structures with information about detected fields.
%   <options>       Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties    Values
%   ----------------------------------------------------------------------------------------------
%    'searchWidth'  If map is not perfect, but contains NaN values along borders, then
%                   search for border pixels can have NaNs. To mitigate this, we check
%                   searchWidth rows/columns near border and if the closest to the border pixel
%                   equals to NaN, we search for first non-NaN value in searchWidth rows-columns.
%                   This argument is optional and default value is 8.
%    'walls'        Definition of walls along which the border score is calculated. Provided by
%                   a string which contains characters that stand for walls:
%                   T - top wall (we assume the bird-eye view on the arena)
%                   R - right wall
%                   B - bottom wall
%                   L - left wall
%                   Characters are case insensitive. Default value is 'TRBL' meaning that border
%                   score is calculated along all walls. Any combination is possible, e.g.
%                   'R' to calculate along right wall, 'BL' to calculate along two walls, e.t.c.
%   ==============================================================================================
%   coverage        Border coverage, ranges from 0 to 1.
%
function coverage = borderCoverage(fields, varargin)
    coverage = 0;

    inp = inputParser;
    defaultSearchWidth = 8;
    defaultWalls = 'TRBL'; % top, right, bottom, left

    checkSearchWidth = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    checkWalls = @internCheckWalls;

    addRequired(inp, 'fields');
    addParameter(inp, 'searchWidth', defaultSearchWidth, checkSearchWidth);
    addParameter(inp, 'walls', defaultWalls, checkWalls);

    inp.KeepUnmatched = true;
    parse(inp, fields, varargin{:});

    % get parsed results
    walls = inp.Results.walls;
    searchWidth = inp.Results.searchWidth;

    % find out what walls are present
    dict = {'B', 'L', 'R', 'T'}; % already sorted
    u = unique(walls);
    t = ['^' sprintf('%c{0,%d}', [u; histc(walls, u)]) '$'];
    wallsIdx = find(~cellfun('isempty', regexpi(dict, t))); % indices in dict of present walls

    for i=1:length(fields)
        curField = fields(i).map;

        % coverage for left wall
        if any(ismember(wallsIdx, 2))
            aux_map = curField(:, 1:searchWidth);
            [covered, norm] = wall_field(aux_map);
            if (covered/norm > coverage)
                coverage = covered/norm;
            end
        end

        % coverage for right wall
        if any(ismember(wallsIdx, 3))
            aux_map = curField(:, end:-1:end+1-searchWidth);
            [covered, norm] = wall_field(aux_map);
            if (covered/norm > coverage)
                coverage=covered/norm;
            end
        end

        % coverage for bottom wall, since we are dealing with data that came from a camera
        % 'bottom' is actually at the top of the matrix curField.
        if (any(ismember(wallsIdx, 1)))
            aux_map = curField(1:searchWidth, :)';
            [covered, norm] = wall_field(aux_map);
            if (covered/norm > coverage)
                coverage=covered/norm;
            end
        end

        % coverage for top wall
        if (any(ismember(wallsIdx, 4)))
            aux_map = curField(end:-1:end+1-searchWidth, :)';
            [covered, norm] = wall_field(aux_map);
            if (covered/norm > coverage)
                coverage = covered/norm;
            end
        end
    end
end

% 'covered' pixels will have distance to border 0.
% Essentially we need to calculate number of elements,
% that equal to zero and take NaNs into account.
function [covered, norm] = wall_field(map)
    ly = size(map, 1);

    D = bwdist(map);
    nanIndices = find(isnan(map(:, 1)));
    numNans = length(nanIndices);
    for i = 1:numNans
        testRow = nanIndices(i);
        nonNan = find(~isnan(map(testRow, :)), 1, 'first');
        if ~isempty(nonNan)
            numNans = numNans - 1;
            D(testRow, 1) = D(testRow, nonNan);
        end
    end

    norm = ly - numNans;
    covered = nansum(D(:, 1) == 0) - numNans;
end

% check input argument that defines walls
% We need to figure out if argument contains one of the characters
% from the dictionary.
% Rise descriptive exception in case of an error.
function res = internCheckWalls(walls)
    res = true;
    if ~ischar(walls)
        error('BNT:args:notChar', 'Argument is not a string.')
    end
    if length(walls) > 4
        error('BNT:args:length', 'Argument can not exceed 4 characters (default ''TLBR'')');
    end

    dict = {'T', 'R', 'B', 'L'};

    % http://stackoverflow.com/questions/19343339/matlab-find-the-indices-of-a-cell-array-of-strings-with-characters-all-containe
    u = unique(walls);
    t = ['^' sprintf('%c{0,%d}', [u; histc(walls, u)]) '$'];
    s = cellfun(@sort, dict, 'uni', 0);
    idx = find(~cellfun('isempty', regexp(s, t)), 1);
    %res = ~isempty(idx); % if it is not empty, then we have characters from dict in walls
    if isempty(idx)
        error('BNT:args:noValidChars', 'There is no information about the walls in the argument.');
    end
end