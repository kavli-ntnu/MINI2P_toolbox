% Locate place fields in a firing map.
%
% Identifies place fields in 2D firing map. Placefields are identified by using
% an adaptive threshold. The idea is that we start with a peak value as the
% threshold. Then we gradually decrease the threshold until the field area doesn't
% change any more or the area explodes (this means the threshold is too low).
%
%  USAGE
%   [fieldsMap, fields] = analyses.palcefieldAdaptive(map, <options>)
%   map           Firing rate map either structure obtained using <a href="matlab:help analyses.map">analyses.map</a>
%                 or a matrix that represents firing map.
%   <options>     optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties     Values
%   -------------------------------------------------------------------------
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
%    'peakCoords'   Matrix Nx2, [x y] or [col row]. If provided, then the function will
%                   not search for peaks, but will use provided values as
%                   peak coordinates. Might be useful if you are interested
%                   in locating fields in some predefined region. Default
%                   value is [].
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
%       [~, fields] = analyses.placefield(map, 'minBins', 5, 'minPeak', 0.1, 'pos', pos);
%       fieldPos = pos(fields(1).posInd, :);
%
function [fieldsMap, fields] = placefieldAdaptive(map, varargin)
    inp = inputParser;
    defaultMinPeak = 1.;
    defaultMinBins = 9;
    defaultBinWidth = 1;
    defaultPos = [];
    defaultMinMeanRate = 0.;
    defaultPeakCoords = [];

    % input argument check functions
    checkScalarZero = @(x) helpers.isiscalar(x, '>=0');
    checkDScalar = @(x) helpers.isdscalar(x, '>0');
    checkDScalarZero = @(x) helpers.isdscalar(x, '>=0');

    % fill input parser object
    addRequired(inp, 'map');
    addParameter(inp, 'minBins', defaultMinBins, checkScalarZero);
    addParameter(inp, 'minPeak', defaultMinPeak, checkDScalarZero);
    addParameter(inp, 'binWidth', defaultBinWidth, checkDScalar);
    addParameter(inp, 'minMeanRate', defaultMinMeanRate, checkScalarZero);
    addParameter(inp, 'pos', defaultPos, @(x) ismatrix(x) && size(x, 2) >= 3);
    addParameter(inp, 'peakCoords', defaultPeakCoords, @(x) ismatrix(x) && ...
        size(x, 2) >= 2 && size(x, 1) >= 1);

    parse(inp, map, varargin{:});

    % get parsed arguments
    minPeak = inp.Results.minPeak;
    minBins = inp.Results.minBins;
    binWidth = inp.Results.binWidth;
    pos = inp.Results.pos;
    minMeanRate = inp.Results.minMeanRate;
    peakCoords = inp.Results.peakCoords;

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

    % Define structuring element disc
    diskSize = 1;
    se = strel('disk', diskSize);

    Ie = imerode(map, se);
    Iobr = imreconstruct(Ie,  map);

    if isempty(peakCoords)
        % this gives regional maxima that can be represented by several pixels,
        % we need maxima consisting of only a single pixel
        regionalMaxMap = imregionalmax(Iobr);

        stats = regionprops(logical(regionalMaxMap), 'Centroid');
        peaks = round(cat(1, stats.Centroid));
    else
        peaks = peakCoords;
    end

    % regionprops returns centroid values in [x y] coordiantes, whereas
    % regular matrix indices are [y x].
    colAxis = 1;
    rowAxis = 2;

    % ensure that there are no peaks off the map (due to rounding)
    peaks(peaks < 1) = 1; % [x, y]
    peaks(peaks(:, colAxis) > size(map, 2), colAxis) = size(map, 2);
    peaks(peaks(:, rowAxis) > size(map, 1), rowAxis) = size(map, 1);

    peakLinInd = sub2ind(size(map), peaks(:, rowAxis), peaks(:, colAxis));

    % obtain peak values. We get them from map (instead of Iobr) in order
    % to get the real values. It may result in inconsistency.
    foundPeaks = map(peakLinInd);

    % remove peaks that have smaller rate than minPeak
    selected = foundPeaks < minPeak;
    peakLinInd(selected) = [];
    peaks(selected, :) = [];
    % Counter for the number of fields
    nFields = 0;

    I = Iobr;
    maxValue = max(I(:));
    % exclude peaks with small values from being detected
    I(I < minPeak) = maxValue * 1.5;

    for i = 1:length(peakLinInd)
        otherFields = peakLinInd;
        otherFields(i) = [];
        if isempty(otherFields)
            otherFields = 0;
        end

        % the values are based on trial and error
        usedTh = 0.96;
        [initialChange, ~, area2] = areaChange(I, peakLinInd(i), 0.96, 0.94, otherFields);
        if isempty(initialChange)
            for j = 0.97:0.01:1
                [initialChange, area1, area2, firstPixels] = areaChange(I, peakLinInd(i), j, j-0.01, otherFields);
                if ~isempty(initialChange)
                    usedTh = j-0.01;
                    break;
                end
            end
            if isempty(initialChange) && ~isempty(area1)
                I(firstPixels) = maxValue * 1.5;
                fieldsMap(firstPixels) = 1;
                continue;
            end
            if isempty(initialChange)
                % failed to extract the field
                continue;
            end
        end
        pixelList = expandField(I, peakLinInd(i), initialChange, area2, otherFields, usedTh);
        if isempty(pixelList)
            [~, pixelList] = areaForThresholdLabel(I, peakLinInd(i), usedTh+0.01, otherFields);
        end

        I(pixelList) = maxValue * 1.5;
        fieldsMap(pixelList) = 1;
    end

    cc = bwconncomp(fieldsMap, 4);
    stats = regionprops(cc, 'Area', 'BoundingBox', 'Centroid', 'PixelIdxList');
    for i = 1:length(stats)
        linInd = stats(i).PixelIdxList;
        [fieldPeak, peakLinearInd] = max(map(linInd));
        [r, c] = ind2sub(size(map), linInd);
        meanRate = nanmean(map(linInd));
        [pr, pc] = ind2sub(size(map), linInd(peakLinearInd));

        if length(r) <= minBins || fieldPeak < minPeak || meanRate < minMeanRate
            continue;
        end

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

function pixelList = expandField(I, peakLinInd, initialChange, initialArea, otherFields, initialTh)
    pixelList = [];
    lastChange = initialChange;
    lastArea = initialArea;
    lastPixels = [];
    numNotChanging = 0;
    numDecrease = 0;

    for i = initialTh:-0.02:0.2
        [area, pixels, isBad] = areaForThresholdLabel(I, peakLinInd, i, otherFields);
        if isempty(area) || isBad
            pixelList = lastPixels;
            break;
        end

        curChange = area / lastArea * 100;
        if curChange < 100
            numDecrease = numDecrease + 1;
        else
            numDecrease = 0;
        end

        if floor(curChange / initialChange) <= 2 && numDecrease < 3
            if curChange == 100
                numNotChanging = numNotChanging + 1;
            else
                numNotChanging = 0;
            end

            if numNotChanging < 10
                lastChange = curChange;
                lastArea = area;
                lastPixels = pixels;
                continue;
            end
        end

        % first if our threshold
        pixelList = lastPixels;
        break;
    end
    if any(lastPixels == peakLinInd)
        % good field
        pixelList = lastPixels;
    end
end

function [ar, pixels, isBad] = areaForThresholdLabel(I, curLinInd, th, otherFields)
    % algorithm based on bwlabeln
    ar = [];
    isBad = false;

    peakValue = I(curLinInd);
    thValue = peakValue * th;
    mask = I >= thValue;

    L = bwlabel(mask, 4);
    pixels = find(L == L(curLinInd));
    % this is needed for Euler number calculation
    bw = padarray(L == L(curLinInd), [1 1]);

    lut = 4*[0 0.25 0.25 0 0.25 0 -.5 -0.25 0.25 -0.5  0 -0.25 0 -0.25 -0.25 0] + 2;
    weights = bwlookup(bw, lut);
    % not exactly the eulerNumber, but enough for us
    eulerNumber = sum(weights(:), 1, 'double') - 2*numel(bw);
    if eulerNumber <= 0
        pixels = [];
        return;
    end

    ar = numel(pixels);
    %isBad = ~isempty(intersect(pixels, otherFields));
    rr = pixels - otherFields';
    isBad = ~all(rr(:));
end

function [ac, area1, area2, firstPixels, secondPixels] = areaChange(I, curLinInd, first, second, otherFields)
    ac = [];
    secondPixels = [];
    area2 = [];

    [area1, firstPixels, isBad] = areaForThresholdLabel(I, curLinInd, first, otherFields);
    if isempty(area1) || isBad
        return;
    end
    [area2, secondPixels, isBad] = areaForThresholdLabel(I, curLinInd, second, otherFields);
    if isempty(area1) || isBad
        return;
    end

    ac = area2 / area1 * 100;
end
