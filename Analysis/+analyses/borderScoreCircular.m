% Calculate border score for a firing map recorded in a circular arena
%
% Calculates a border score for a firing rate map according to the article "Representation of Geometric
% Borders in the Entorhinal Cortex" by Solstad et. al. (Science 2008).
% Border score ranges from -1 to +1 with +1 being "a perfect border cell". If the firing map contains
% no firing fields, then the returned score is -1.
% The score reflects not only how close a field is to a border and how big the coverage of this field is,
% but also it reflects spreadness of a field. The best border score (+1) will be calculated for a thin
% line (1 px, bin) that lies along the wall and fully covers this wall.
% There are some differences between calculating border score for rectangular and circular environment:
% 1. walls are defined differently. There are no obvious top, bottom, left, right walls. User should define
%    walls she is interested in.
% 2. Normalization by mean firing distance is done through a sigmoid function of coverage. This means
%    that mean firing distance is weighted by sigmoid(coverage).
%
%  USAGE
%   score = analyses.borderScoreCircular(map, fieldsMap, <options>)
%   map         Matrix, a 2D firing map in polar coordinates.  This
%               means that x-axis should contain values in range 1-360, y-axis range
%               is 1 to radius of the circular arena.
%   fieldsMap   Matrix of the same size as map with only non-zero elements
%               belonging to detected fields. This is an output from
%               function analyses.placefield.
%   <options>   Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties    Values
%   ----------------------------------------------------------------------------------------------
%    'searchWidth'  If map is not perfect, but contains NaN values along borders, then
%                   search for border pixels can have NaNs. To mitigate this, we check
%                   searchWidth rows/columns near border and if the closest to the border pixel
%                   equals to NaN, we search for first non-NaN value in searchWidth rows/columns.
%                   This argument is optional and default value is 8.
%    'walls'        Definition of walls along which the coverage is calculated. Provided by
%                   a 2D matrix Nx2. Each row defines 1 wall. Each row should have two values in degrees:
%                   1. angle at which wall starts.
%                   2. angle at which wall ends.
%                   Default value is [1 360], which represents a single wall around the whole circular
%                   environment. Note that you can not use 0 as an angle.
%   ==============================================================================================
%   score       Border score. Ranges from -1 to +1 (if there are no fields, then the score is -1).
%
%  SEE
%   See analyses.placefield
%
function score = borderScoreCircular(map, fieldsMap, varargin)
    score = -1;
    if ~any(fieldsMap >= 1)
        % no fields
        return;
    end

    coverage = analyses.borderCoverageCircular(fieldsMap, varargin{:});

    fieldsFiringMap = fieldsMap;
    fieldsFiringMap(fieldsFiringMap > 1) = 1; % different fields are assigned different numbers,
                                              % make them all 1.
    fieldsFiringMap = map .* fieldsFiringMap;

    fdist = weighted_firing_distance(fieldsFiringMap);
    % construct a sigmoid activation function
    x = (0:0.01:1) * 10;
    y = sigmf(x, [-2 4]);
    % get fdist activation value through interpolation of coverage
    % the idea is that with higher coverage fdist should be lower.
    fdistC = interp1(x, y, coverage*10);
    fdist = fdist * fdistC;

    score = (coverage - fdist)/(coverage + fdist);
    % fprintf('Coverage %f, d: %f, score = %f\n', coverage, fdist, score);
end

% This function differs from the version for rectangular arenas
% WFD plays an important role. It should not be a big value, but it also should not be a small value.
% Probably, a constant can be used. WFD regulates the score and controls 'the diverseness' of the score.
% If we take one cell and calculate it's score with different walls (so that coverage is different), WFD
% will be constant, but the score will be different. For example, consider a cell that covers range 260..360
% almost perfectly. Score for such wall can be 0.558957. The same cell, but with wall [1 360] yields score
% 0.072770.
function wfd = weighted_firing_distance(map)
    [ly, lx] = size(map);
    [~, my] = meshgrid(1:lx, 1:ly);
    % Since map is in polar coordinates, y-axis values are distances already;
    distance_matrix = flipud(my);

    wfd = nanmean(nanmean(map .* distance_matrix));

    % x-axis is one big wall (again because map is in polar coordinates)
    wfd = wfd / ly;
end
