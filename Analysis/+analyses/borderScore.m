% Calculate border score for a firing map
%
% Calculates a border score for a firing rate map according to the article "Representation of Geometric
% Borders in the Entorhinal Cortex" by Solstad et. al. (Science 2008).
% Border score ranges from -1 to +1 with +1 being "a perfect border cell". If the firing map contains
% no firing fields, then the returned score is -1.
% The score reflects not only how close a field is to a border and how big the coverage of this field is,
% but also it reflects spreadness of a field. The best border score (+1) will be calculated for a thin
% line (1 px, bin) that lies along the wall and fully covers this wall.
%
%  USAGE
%   score = analyses.borderScore(map, fieldsMap, fields, <options>)
%   map         A 2D firing map, not a struct you get from analyses.map, but actual values.
%   fieldsMap   Matrix of the same size as map with only non-zero elements
%               belonging to detected fields. This is an output from
%               function analyses.placefield.
%   fields      Array of structures with information about detected fields.
%               This is an output from function analyses.placefield.
%   <options>   Optional list of property-value pairs (see table below)
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
%   score       Border score. Ranges from -1 to +1 (if there are no fields,
%               then score is -1).
%
%  SEE
%   See analyses.placefield
%
function score = borderScore(map, fieldsMap, fields, varargin)
    score = -1;

    if ~isempty(fields)
        coverage = analyses.borderCoverage(fields, varargin{:});

        fieldsFiringMap = fieldsMap;
        fieldsFiringMap(fieldsFiringMap > 0) = 1; % different fields are assigned different numbers,
                                                  % make them all 1.
        fieldsFiringMap = map .* fieldsFiringMap;

        fdist = weighted_firing_distance(fieldsFiringMap);
        score = (coverage - fdist)/(coverage + fdist);
    end
end

function wfd = weighted_firing_distance(map)
    % normalized firing map. Normalization is done by
    % the sum of firing rates of all pixels belonging to all fields
    map = map/nansum(map(:));

    [ly, lx] = size(map);
    [mx, my] = meshgrid(1:lx, 1:ly);

    % alternative code for distance matrix is:
    % map = zeros(ly, lx);
    % map(1, :) = 1;
    % map(end, :) = 1;
    % map(:, 1) = 1;
    % map(:, end) = 1;
    % D = bwdist(map);
    % but this code is slower.

    distance_matrix = min(min(my,mx), min(flipud(my), fliplr(mx)));

    wfd = nansum(nansum(map .* distance_matrix));

    % normalization by half of the smallest arena size min(ly, lx)/2.
    wfd = (2 * wfd) / min(ly, lx);
end
