% Calculate gridness score for an autocorrelogram
%
% Calculates a gridness score by expanding a circle around the centre field and
% calculating a correlation value of the expanded circle with it's rotated versions.
% The expansion is done up until the smallest side of the autocorrelogram.
% Can also calculate grid statistics.
%
% Gridness score value by itslef is calculated as a maximum over sliding mean
% of expanded circles. The widtg of the sliding window is given by variable
% numGridnessRadii. This is done in order to keep the values the same with
% historical development of gridness score.
%
% NB! The function assumes that the provided input is indeed an autocorrelgoram.
% This means that the center pixel must have the maximum value. This must be true
% for both square and rectangular autocorrelogram. This will be important if you
% manually create square matrix from a rectangular one. Imagine you have a rectangular
% matrix aCorr of size 309x305, aCorr(155, 153) == 1. If you make it rectangular
% like  this `aSquare = aCorr(1:305, :);`, then the centre pixel will not have
% the highest values, i.e. `aSquare(153, 153) ~= 1`. You should do it like this
% `aSquare = aCorr(3:309-2, :)`, then `aSquare(153, 153)` yields `1`.
%
%  USAGE
%   [score, <stats>] = analyses.gridnessScore(aCorr, <options>)
%   aCorr       A 2D autocorrelogram. aCorr can be square or rectangular. See the note (NB) above!
%   <options>   Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties         Values
%   ----------------------------------------------------------------------------------------------
%    'minOrientation'   Value of minimal difference of inner fields orientation (in degrees). If
%                       there are fields that differ in orientation for less than minOrientation,
%                       then only the closest to the centre field are left. Default value is 15.
%   'debug'             True or False. If set to True, the function produces some debug output.
%   ==============================================================================================
%   score       Gridness score. Ranges from -2 to 2. 2 is more a theoretical bound for a perfect grid.
%               More practical value is around 1.3.
%   stats       If this variable is requested, then it is a structure with the following statistics:
%       spacing         3-element vector with distances from the centre field to neighbour fields.
%       orientation     3-element vector with orientations between the centre field and neighbour fields.
%       ellipse         Ellipse fitted to the grid. Contains the centre, radii and orientation in
%                       radians, stored as [Cx, Cy, Rx, Ry, theta].
%       ellipseTheta    Radius of the ellipse in degrees wrapped in range [0..180].
%
function [gscore, varargout] = gridnessScore(aCorr, varargin)
    nout = max(nargout, 1) - 1;
    inp = inputParser;
    defaultMinOrientation = 30;
    gridStat.orientation = [];
    gridStat.spacing = [];
    gridStat.ellipse = [];
    gridStat.ellipseTheta = nan;

    % input argument check functions
    checkDScalar = @(x) helpers.isdscalar(x, '>0');

    % fill input parser object
    addRequired(inp, 'aCorr');
    addParameter(inp, 'minOrientation', defaultMinOrientation, checkDScalar);
    addParameter(inp, 'debug', false);

    parse(inp, aCorr, varargin{:});

    % get parsed arguments
    minOrientation = inp.Results.minOrientation;
    isDebug = inp.Results.debug;

    halfSize = ceil(size(aCorr)/2);
    half_height = halfSize(1);
    half_width = halfSize(2);
    aCorrRad = min(halfSize);
    aCorrSize = size(aCorr);

    if aCorrSize(1) == 1 || aCorrSize(2) == 1
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
        end
        return;
    end

    % contourc is efficient if aCorr is normalized
    maxValue = max(max(aCorr));
    if maxValue ~= 1
        aCorr = aCorr / maxValue;
    end

    cFieldRadius = findCentreRadius(aCorr, half_width, half_height);
    if isDebug
        fprintf('Center radius is %f\n', cFieldRadius);
    end
    if cFieldRadius == 0 || cFieldRadius == 1 || cFieldRadius == -1 || cFieldRadius >=min(halfSize)
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
        end
        return;
    end

    % Meshgrid for expanding circle
    [rr, cc] = meshgrid(1:size(aCorr, 2), 1:size(aCorr, 1));

    % Define iteration radius step size for the gridness score
%     if cFieldRadius>=aCorrRad  % modified by Weijian Zong,20201012
%     cFieldRadius=aCorrRad;  % modified by Weijian Zong,20201012
%     else% modified by Weijian Zong,20201012
%     end % modified by Weijian Zong,20201012
    radSteps = cFieldRadius:aCorrRad;
    radSteps(1) = [];
    numSteps = length(radSteps);

    GNS = zeros(numSteps, 2);
    rotCorr = zeros(1, 5);
    rotAngles_deg = 30*(1:5);

    % aCorr is rotated outside the loop for speed
    rotatedCorr = cell(1, length(rotAngles_deg));
    for i = 1:length(rotAngles_deg)
        rotatedCorr{i} = imrotate(aCorr, rotAngles_deg(i), 'bilinear', 'crop');
    end

    mainCircle = sqrt((cc - half_height).^2 + (rr - half_width).^2);
    innerCircle = mainCircle > cFieldRadius;

    % Define expanding ring of autocorrellogram and do x30 correlations
    for i = 1:numSteps
        ind = (innerCircle & (mainCircle < radSteps(i)));
        tempCorr = reshape(aCorr(ind), 1, [])';
        for j = 1:5
            rotatedCircle = reshape(rotatedCorr{j}(ind), 1, [])';
            rotCorr(j) =  corr(tempCorr, rotatedCircle);
            if isDebug
                fprintf('Step %u, angle %u, corr value %f\n', i-1, rotAngles_deg(j), rotCorr(j));
            end
        end
        GNS(i, 1) = min(rotCorr([2, 4])) - max(rotCorr([1, 3, 5]));
        GNS(i, 2) = radSteps(i);
    end

    % Find the biggest gridness score and radius
    [~, gscoreLoc] = max(GNS(:, 1));
    % See function help about numGridnessRadii
    numGridnessRadii = 3;
    numStep = numSteps - numGridnessRadii;
    if numStep < 1
        numStep = 1;
    end

    if numStep == 1
        gscore = nanmean(GNS(:, 1));
    else
        meanGridnessArray = zeros(numStep, 1);
        for ii = 1:numStep
            meanGridnessArray(ii) = nanmean(GNS(ii:ii + numGridnessRadii-1, 1));
        end

        [gscore, gInd] = max(meanGridnessArray);
        gscoreLoc = gInd + (numGridnessRadii-1)/2;
    end

    varargout{4} = radSteps(gscoreLoc);

    % Return if we do not need to calculate grid statistics
    if nout < 1
        return;
    end

    %% Calculate gridness score statistics
    bestCorr = (mainCircle < radSteps(gscoreLoc) * 1.25) .* aCorr;
    regionalMaxMap = imregionalmax(bestCorr, 4);
    se = strel('square', 3);
    im2 = imdilate(regionalMaxMap, se); % dilate map to eliminate fragmentation
    cc = bwconncomp(im2, 8);
    stats = regionprops(cc, 'Centroid');

    if length(stats) < 5
        warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

        varargout{1} = gridStat;
        varargout{2} = cFieldRadius;
        varargout{3} = nan(6, 2);
        varargout{4} = radSteps(gscoreLoc);
        return;
    end

    allCoords = [stats(:).Centroid];
    centresOfMass(:, 1) = allCoords(1:2:end);
    centresOfMass(:, 2) = allCoords(2:2:end);

    % Calculate orientation for each field relative to the centre field
    orientation = (atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)
    peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    zeroInd = find(orientation == 0, 1);
    orientation(zeroInd) = []; % remove zero value, so that we do not have a side effect with minOrientation
    stats(zeroInd) = [];
    peaksToCentre(zeroInd) = [];
    centresOfMass(zeroInd, :) = [];

    % filter fields that have similar orientation
    orientDistSq = circ_dist2(orientation);
    closeFields = abs(orientDistSq) < deg2rad(minOrientation);
    [rows, cols] = size(closeFields);
    closeFields(1:(rows+1):rows*cols) = 0; % assign zero to diagonal elements
    closeFields(tril(true(rows))) = 0; % assign zero to lower triangular of a matrix. Matrix is
                                       % symmetric and we do not need these values.
    [rows, cols] = find(closeFields); % find non-empty elements, they correspond to indices of close fields
    if ~isempty(rows)
        indToDelete = zeros(1, length(rows));
        for i = 1:length(rows)
            % fieldPeaks = [fields([rows(i) cols(i)]).peakX; fields([rows(i) cols(i)]).peakY];
            % peaksToCentre = sqDistance(fieldPeaks, [half_width; half_height]);
            if peaksToCentre(rows(i)) > peaksToCentre(cols(i))
                indToDelete(i) = rows(i);
            else
                indToDelete(i) = cols(i);
            end
        end
        indToDelete = unique(indToDelete);
        stats(indToDelete) = [];
        peaksToCentre(indToDelete) = [];

        if length(stats) < 4
            warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

            varargout{1} = gridStat;
            varargout{2} = cFieldRadius;
            varargout{3} = nan(6, 2);
            varargout{4} = radSteps(gscoreLoc);
            return;
        end

        allCoords = [stats(:).Centroid];
        clear centresOfMass;
        centresOfMass(:, 1) = allCoords(1:2:end);
        centresOfMass(:, 2) = allCoords(2:2:end);
    end

    % % get fields peak coordinates
    % fieldPeaks = zeros(length(stats), 2);
    % for i = 1:length(stats)
    %     [~, maxInd] = max(bestCorr(stats(i).PixelIdxList));
    %     fieldPeaks(i, :) = stats(i).PixelList(maxInd, :);
    % end
    % % fieldPeaks = [fields(:).peakX; fields(:).peakY]'; %
%     peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    [~, sortInd] = sort(peaksToCentre);
    stats = stats(sortInd);
    centresOfMass = centresOfMass(sortInd, :);

    % leave only 6 closest neighbours (if available)
    if length(stats) > 5
%         stats = stats(1:6);
        centresOfMass = centresOfMass(1:6, :);
    else
%         stats = stats(1:end);
        centresOfMass = centresOfMass(1:end, :);
    end

    % centresOfMass = [fields(:).x; fields(:).y]';
    % Calculate orientation for each field relative to the centre field
    orientation = rad2deg(atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)

    % Calculate distances between centre of masses for each field and the centre field
    spacing = sqrt((centresOfMass(:, 1) - half_width).^2 + (centresOfMass(:, 2) - half_height).^2);

%     % Plot grid polygon points
%     figure, plot.colorMap(bestCorr), hold on;
%     plot(centresOfMass(:, 1), centresOfMass(:, 2), '+k', 'markersize', 8);

    ell = general.fitEllipse(centresOfMass(:, 1), centresOfMass(:, 2));
    ellipseTheta = rad2deg(general.wrap(ell(end)) + pi);
%     drawEllipse(ell, 'linewidth', 2, 'color', [1 1 1]);

    % Determine axes orientation, spacing and deviation
    [~, bBC] = sort(abs(orientation));
    [~, bBC2] = sort(abs(orientation - orientation(bBC(1))));

    % leave only three values, because autocorrelogram is symmetric
    orientation = orientation(bBC2(1:3));
    spacing = spacing(bBC2(1:3));
    [orientation, orientSortInd] = sort(orientation);

    spacing = spacing(orientSortInd);

    gridStat.orientation = orientation;
    gridStat.spacing = spacing;
    gridStat.ellipse = ell;
    gridStat.ellipseTheta = ellipseTheta;

    varargout{1} = gridStat;
    varargout{2} = cFieldRadius;
    varargout{3} = centresOfMass;
    varargout{4} = radSteps(gscoreLoc);
end

function D = sqDistance(X, Y)
    D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
end

function cFieldRadius = findCentreRadius(aCorr, half_width, half_height)
    cFieldRadius = 0;

    % Search for fields only around the centre
    peakCoords = [half_width, half_height];
    [~, fields] = analyses.placefieldAdaptive(aCorr, 'minPeak', 0, 'minBins', 2, ...
        'peakCoords', peakCoords);
    if isempty(fields)
        return;
    end

    peakLoc = [fields(:).peakX; fields(:).peakY]';

    % get all distances and check two minimums of them
    allDistances = sqDistance(peakLoc', [half_width half_height]'); % point should be in format [x y]
    [~, sortIndices] = sort(allDistances);
    if length(sortIndices) >= 2
        % this is a bit leagacy code. By using peakCoords with placefieldAdaptive, we
        % should always get just a single field. Keeping it as I did not test the version without it properly.
        twoMinIndices = sortIndices(1:2);

        if abs(allDistances(twoMinIndices(1)) - allDistances(twoMinIndices(2))) < 2
            % two fields with close middle points. Let's select one with minimum square
            [~, minInd] = min([fields(twoMinIndices).area]);
            closestFieldInd = twoMinIndices(minInd);
        else
            closestFieldInd = twoMinIndices(1); % get the first minimum
        end
    else
        closestFieldInd = sortIndices(1);
    end
    cFieldRadius = floor(sqrt(fields(closestFieldInd).area / pi));
end