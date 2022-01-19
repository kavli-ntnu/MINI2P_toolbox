% Locate place fields in a firing map.
%
% Identifies the place fields in 2D firing map. First map is converted to binary image
% by applying a global threshold. During the second step connected regions are labeled
% by function bwlabel. Final step includes filtering of identified regions and verification
% that they are valid place fields.
%
%  USAGE
%   [fieldsMap, fields] = analyses.palcefield(map, <options>)
%   map           Firing rate map either structure obtained using <a href="matlab:help analyses.map">analyses.map</a>
%                 or a matrix that represents firing map.
%   <options>     optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties     Values
%   -------------------------------------------------------------------------
%    'threshold'    value above threshold*peak belong to a field (default = 0.2).
%    'minBins'      Minimum number of bins in a place field, i.e. total number
%                   of bins of place field or an area of place field. Fields with
%                   fewer bins are not considered as place fields. Remember to
%                   adjust this value when you change the bin width.
%                   (default = 9).
%    'minPeak'      peaks smaller than this value are considered spurious
%                   and ignored (default = 1). Peak is normally a rate, however
%                   it's units not necessary are Hz.
%    'binWidth'     Width of the bins in cm. It is used to calculate field size.
%                   (default = 1).
%    'pos'          Position samples. Used to calculate posInd. If not provided,
%                   then posInd will be an empty matrix.
%    'minMeanRate'  Fields with mean rate smaller that this value are ignored (default = 0).
%   =========================================================================
%
%  OUTPUT
%
%   fieldsMap       Matrix of the same size as firing map. The elements of fieldsMap
%                   are integer values greater or equal to 0. The elements labeled 0
%                   are the background (not a place field). The pixels labeled 1 make up
%                   first field; the pixels labeled 2 make up a second object; and so on.
%
%   fields          Structure with information about each field. Structure fields are:
%       row         vector of rows that constitute this field;
%       col         vector of columns that constitute this field;
%       area        area of field measured by regionprops function, i.e. number of bins in the field;
%       bbox        bounding box of the field;
%       peak        field peak value;
%       size        field size in cm; Calculated with help of binWidth argument.
%       x, y        field centre of mass point;
%       meanRate    mean firing rate;
%       PixelIdxList    linear list of indices that constitute the field.
%       peakX       X-coordinate of field's peak
%       peakY       Y-coordinate of field's peak
%       map         Binary map of the field. It has the same size as firing map. The elements of 'map'
%                   equal to one if they belong to the field and are zeros otherwise.
%       posInd      Indices of position samples that correspond to this field. In case position sample
%                   matrix is of size Nx5 ([t x y x1 y1]), posInd corresponds to the left most position
%                   columns ([x y]). If pos argument is not provided, then posInd will be an empty matrix.
%
%  EXAMPLES
%
%   1. To obtain positions of the first field run:
%       [~, fields] = analyses.placefield(map, 'threshold', 0.3, 'minBins', 5, 'minPeak', 0.1, 'pos', pos);
%       fieldPos = pos(fields(1).posInd, :);
%

% Note to developers:
% In order to improve speed and capture all fields the map is processed in two passes. First
% we extract fields based solely on global peak value. This should be enough for most of the cases.
% Secondly we check if extracted fields maxima is the same as results of imregionalmax. If they are
% different, then we perform another extraction, which is slower.
%
function [fieldsMap, fields] = placefield(map, varargin)
    inp = inputParser;
    defaultMinPeak = 1.;
    defaultThreshold = 0.2;
    defaultMinBins = 9;
    defaultBinWidth = 1;
    defaultPos = [];
    defaultMinMeanRate = 0;

    % input argument check functions
    checkThreshold = @(x) helpers.isdscalar(x, '>=0', '<=1');
    checkScalarZero = @(x) helpers.isiscalar(x, '>=0');
    checkDScalar = @(x) helpers.isdscalar(x, '>0');
    checkDScalarZero = @(x) helpers.isdscalar(x, '>=0');

    % fill input parser object
    addRequired(inp, 'map');
    addParameter(inp, 'threshold', defaultThreshold, checkThreshold);
    addParameter(inp, 'minBins', defaultMinBins, checkScalarZero);
    addParameter(inp, 'minPeak', defaultMinPeak, checkDScalarZero);
    addParameter(inp, 'binWidth', defaultBinWidth, checkDScalar);
    addParameter(inp, 'minMeanRate', defaultMinMeanRate, checkScalarZero);
    addParameter(inp, 'pos', defaultPos, @(x) ismatrix(x) && size(x, 2) >= 3);

    parse(inp, map, varargin{:});

    % get parsed arguments
    minPeak = inp.Results.minPeak;
    threshold = inp.Results.threshold;
    minBins = inp.Results.minBins;
    binWidth = inp.Results.binWidth;
    pos = inp.Results.pos;
    minMeanRate = inp.Results.minMeanRate;

    originalMap = [];
    if isstruct(map)
        originalMap = map;
        map = map.z;
    end

    fieldsMap = zeros(size(map));
    fields = struct('row', {}, 'col', {}, ...
        'size', {}, 'peak', {}, 'peakX', {}, 'peakY', {}, ...
        'area', {}, 'bbox', {}, ...
        'x', {}, 'y', {}, ...
        'meanRate', {}, 'PixelIdxList', {}, ...
        'map', {}, 'posInd', {} ...
        );

    if isempty(map)
        return
    end

    globalPeak = nanmax(nanmax(map));
    if isnan(globalPeak) || globalPeak == 0
        return;
    end

    mapNans = isnan(map);
    map(mapNans) = 0;
    regionalMaxMap = imregionalmax(map, 4); % obtain all local maxima
    testRegionalMx = zeros(size(map)); % map that will contain located fields maxima

    [ir, ic] = find(regionalMaxMap > 0);            % get peak locations
    foundPeaks = map(sub2ind(size(map), ir ,ic));   % obtain peaks value

    % remove peaks that have smaller rate than minPeak
    selected = foundPeaks < minPeak;
    regionalMaxMap(ir(selected), ic(selected)) = 0;

    % Counter for the number of fields
    nFields = 0;

    binMap = ones(size(map));

    binMap(map < globalPeak*threshold) = 0;
    binMap(mapNans) = 0;

    cc = bwconncomp(binMap, 4);
    stats = regionprops(cc, 'Area', 'BoundingBox', 'Centroid');

    for i = 1:cc.NumObjects
        linInd = cc.PixelIdxList{i};
        [fieldPeak, peakLinearInd] = max(map(linInd));
        [r, c] = ind2sub(size(map), linInd);
        meanRate = nanmean(map(linInd));

        [pr, pc] = ind2sub(size(map), linInd(peakLinearInd));
        % mark this field as visited even if we reject if further down
        testRegionalMx(pr, pc) = 1;

        if fieldPeak < minPeak
            continue;
        end
        if meanRate < minMeanRate
            continue;
        end

        if length(r) >= minBins
            nFields = nFields + 1;

            fields(nFields).row = r;
            fields(nFields).col = c;
            fields(nFields).size = length(r) * binWidth^2;
            fields(nFields).peak = fieldPeak;
            fields(nFields).peakX = pc;
            fields(nFields).peakY = pr;
            fields(nFields).area = stats(i).Area;
            fields(nFields).bbox = stats(i).BoundingBox;
            fields(nFields).PixelIdxList = linInd;

            fields(nFields).x = stats(i).Centroid(1);
            fields(nFields).y = stats(i).Centroid(2);

            fields(nFields).meanRate = meanRate;
            fields(nFields).map = zeros(size(map));
            fields(nFields).map(linInd) = 1;
            fields(nFields).map(mapNans) = nan;

            % fields(nFields).Extent = stats(i).Extent;
            % fields(nFields).Orientation = stats(i).Orientation;
            % fields(nFields).Eccentricity = stats(i).Eccentricity;
            % fields(nFields).MinorAxisLength = stats(i).MinorAxisLength;
            % fields(nFields).MajorAxisLength = stats(i).MajorAxisLength;
            % fields(nFields).Perimeter = stats(i).Perimeter;

            if ~isempty(pos)
                intRect = ceil(fields(nFields).bbox); % bounding box with integers
                if isempty(originalMap)
                    % we do not have position distribution from the map, so assume that it's based
                    % on data min/max values
                    nBins = size(map);
                    limitsX = [nanmin(pos(:, bntConstants.PosX)) nanmax(pos(:, bntConstants.PosX))];
                    limitsY = [nanmin(pos(:, bntConstants.PosY)) nanmax(pos(:, bntConstants.PosY))];
                    xSpace = linspace(limitsX(1), limitsX(2), nBins(2));
                    ySpace = linspace(limitsY(1), limitsY(2), nBins(1));
                else
                    xSpace = originalMap.x;
                    ySpace = originalMap.y;
                end

                xMin = xSpace(intRect(1));
                xMax = xSpace(intRect(1) + intRect(3) - 1);

                yMin = ySpace(intRect(2));
                yMax = ySpace(intRect(2) + intRect(4) - 1);

                posIndX = pos(:, bntConstants.PosX) >= xMin & pos(:, bntConstants.PosX) <= xMax;
                posIndY = pos(:, bntConstants.PosY) >= yMin & pos(:, bntConstants.PosY) <= yMax;
                fields(nFields).posInd = find(posIndX & posIndY);
            else
                fields(nFields).posInd = [];
            end

            fieldsMap(linInd) = nFields;
        end
    end

    if ~isequaln(testRegionalMx, regionalMaxMap)
        map(fieldsMap > 0) = 0; % turn off map values for known fields, prevent field duplicates

        % we have some uncounted fields
        leftPeaksMap = regionalMaxMap - testRegionalMx;
        [ir, ic] = find(leftPeaksMap > 0);                    % get peak locations
        foundPeaks = map(sub2ind(size(leftPeaksMap), ir ,ic));   % obtain peaks value
        mapThresholds = foundPeaks * threshold;
        if sum(foundPeaks) == 0
            return;
        end

        finalMap = zeros(size(map));
        for i = 1:length(foundPeaks)
            binMap = zeros(size(map));
            binMap( (map > mapThresholds(i)) & (map <= foundPeaks(i)) ) = 1;
            binMap(mapNans) = 0;

            [binMap, linInd] = bwselect(binMap, ic(i), ir(i), 4); % leave only region that relates to current field
            if isempty(linInd)
                continue;
            end

            % check for minimum number of bins
            [r, ~] = ind2sub(size(map), linInd);
            if length(r) < minBins
                continue;
            end

            stats = regionprops(binMap, 'Centroid', 'EulerNumber'); % find statistics on field candidates
            distToPeak = sqrt((stats.Centroid(1) - ic(i))^2 + (stats.Centroid(2) - ir(i))^2);
    %         distToPeak = pdist2(stats.Centroid, [ic(i) ir(i)]);
            if distToPeak > 4 || stats.EulerNumber < 1 % we want object without holes (euler number)
                continue;
            end
            finalMap(linInd) = 1;
        end

        cc = bwconncomp(finalMap, 4); % somehow regionprops is not working if finalMap is passed directly
        stats = regionprops(cc, 'Area', 'BoundingBox', 'Centroid');

        for fieldInd = 1:length(stats)
            linInd = cc.PixelIdxList{fieldInd};
            [r, c] = ind2sub(size(map), linInd);
            meanRate = nanmean(map(linInd));
            [peakRate, peakInd] = nanmax(map(linInd));
            [pr, pc] = ind2sub(size(map), linInd(peakInd));

            if peakRate < minPeak
                continue;
            end
            if meanRate < minMeanRate
                continue;
            end

            nFields = nFields + 1;

            fields(nFields).row = r;
            fields(nFields).col = c;
            fields(nFields).size = length(r) * binWidth^2;
            fields(nFields).peak = peakRate;
            fields(nFields).peakX = pc;
            fields(nFields).peakY = pr;
            fields(nFields).area = stats(fieldInd).Area;
            fields(nFields).bbox = stats(fieldInd).BoundingBox;
            fields(nFields).PixelIdxList = linInd;

            fields(nFields).x = stats(fieldInd).Centroid(1);
            fields(nFields).y = stats(fieldInd).Centroid(2);

            fields(nFields).meanRate = meanRate;
            fields(nFields).map = zeros(size(map));
            fields(nFields).map(linInd) = 1;
            fields(nFields).map(mapNans) = nan;

            if ~isempty(pos)
                intRect = ceil(fields(nFields).bbox); % bounding box with integers
                if isempty(originalMap)
                    % we do not have position distribution from the map, so assume that it's based
                    % on data min/max values
                    nBins = size(map);
                    limitsX = [nanmin(pos(:, bntConstants.PosX)) nanmax(pos(:, bntConstants.PosX))];
                    limitsY = [nanmin(pos(:, bntConstants.PosY)) nanmax(pos(:, bntConstants.PosY))];
                    xSpace = linspace(limitsX(1), limitsX(2), nBins(2));
                    ySpace = linspace(limitsY(1), limitsY(2), nBins(1));
                else
                    xSpace = originalMap.x;
                    ySpace = originalMap.y;
                end

                xMin = xSpace(intRect(1));
                xMax = xSpace(intRect(1) + intRect(3) - 1);

                yMin = ySpace(intRect(2));
                yMax = ySpace(intRect(2) + intRect(4) - 1);

                posIndX = pos(:, bntConstants.PosX) >= xMin & pos(:, bntConstants.PosX) <= xMax;
                posIndY = pos(:, bntConstants.PosY) >= yMin & pos(:, bntConstants.PosY) <= yMax;
                fields(nFields).posInd = find(posIndX & posIndY);
            else
                fields(nFields).posInd = [];
            end

            fieldsMap(linInd) = nFields;
        end
    end
end
