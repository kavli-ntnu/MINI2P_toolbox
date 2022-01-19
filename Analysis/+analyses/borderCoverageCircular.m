% Calculate border coverage for detected fields in a circular arena
%
% This function calculates firing map border coverage that is further used
% in calculation of a border score. This function must be used with recordings
% done in a circular environment. See USAGE for details.
%
%  USAGE
%   coverage = analyses.borderCoverageCircular(fieldsMap, <options>)
%   fieldsMap   2D binary matrix that repressents firing properties of the field.
%               fieldsMap must be a polar version of a regular firing rate map. This
%               means that x-axis should contain values in range 1-360, y-axis range
%               is 1 to radius of the circular arena.
%               fieldsMap should contain NaN for unvisited bins. Other possible values:
%               0 for zero firing rate, i.e. area that doesn't belong to a fields.
%               Any positive number marks an area of a field. For example, if there is
%               just a single field, then fieldsMap could consist of 3 values: NaNs, zeros,
%               and ones.
%   <options>   Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties    Values
%   ----------------------------------------------------------------------------------------------
%    'searchWidth'  If map is not perfect, but contains NaN values along borders, then
%                   search for border pixels can have NaNs. To mitigate this, we check
%                   searchWidth rows/columns near border and if the closest to the border pixel
%                   equals to NaN, we search for first non-NaN value in searchWidth rows/columns.
%                   This argument is optional and default value is 8 bins.
%    'walls'        Definition of walls along which the coverage is calculated. Provided by
%                   a 2D matrix Nx2. Each row defines 1 wall. Each row should have two values in degrees:
%                   1. angle at which wall starts.
%                   2. angle at which wall ends.
%                   Default value is [1 360], which represents a single wall around the whole circular
%                   environment. Note that you can not use 0 as an angle.
%   ==============================================================================================
%   coverage    Border coverage, ranges from 0 to 1.
%
function coverage = borderCoverageCircular(fieldsMap, varargin)
    coverage = 0;

    inp = inputParser;
    defaultSearchWidth = 8;
    defaultWalls = [1 360];

    checkSearchWidth = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    checkWalls = @(x) size(x, 2) == 2 && all(x(:) > 0) && all(x(:) <= 360);

    addRequired(inp, 'fieldsMap', @(x) isnumeric(x) && size(x, 2) == 360);
    addParameter(inp, 'searchWidth', defaultSearchWidth, checkSearchWidth);
    addParameter(inp, 'walls', defaultWalls, checkWalls);

    parse(inp, fieldsMap, varargin{:});

    walls = inp.Results.walls;
    searchWidth = inp.Results.searchWidth;
    fieldsMap(fieldsMap > 1) = 1; % make it binary

    for i = 1:size(walls, 1)
        wall = walls(i, :);

        aux_map = fieldsMap(end-searchWidth+1:end, wall(1):wall(2));
        [covered, norm] = wall_field(aux_map);
        coverage = max([covered/norm, coverage]);
    end
end

% 'covered' pixels will have distance to border 0.
% Essentially we need to calculate number of elements,
% that equal to zero and take NaNs into account.
function [covered, norm] = wall_field(map)
    lx = size(map, 2);
    map = flipud(map);

    D = bwdist(map);
    nanIndices = find(isnan(map(1, :)));
    numNans = length(nanIndices);
    for i = 1:numNans
        testColumn = nanIndices(i);
        nonNan = find(~isnan(map(:, testColumn)), 1, 'first');
        if ~isempty(nonNan)
            numNans = numNans - 1;
            D(1, testColumn) = D(nonNan, testColumn);
        end
    end

    norm = lx - numNans;
    covered = nansum(D(1, :) == 0) - numNans;
end
