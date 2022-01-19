% Calculate arena coverage (the amount of space an animal has covered during experiment)
%
% Calculates arena coverage.
%
%  USAGE
%   coverage = analyses.arenaCoverage(pos, binWidth, shape, dimensions)
%   pos         Position samples, matrix of size at least Nx2. Format is either [t x] or [t x y].
%   binWidth    Width of horizontal and vertical bins in cm. If only one value is provided, then
%               the same bin width is used in both directions.
%   shape       Arena shape, integer. One of the values of bntConstants.ArenaShape.
%   dimensions  Vector of arena dimensions, the actual number of elements depends on arena shape.
%
%   coverage    Arena coverage, float in range [0..100].
%
function coverage =  arenaCoverage(pos, binWidth, shape, dimensions)
    inp = inputParser;

    checkPosDimensions = @(x) size(x, 2) >= 2;
    checkBinWidth = @(x) helpers.isdvector(x, '>=0') && length(x) <= 2;
    checkShape = @(x) ismember(x, helpers.ArenaShape.allShapes());
    checkDims = @(x) helpers.isdvector(x, '>=0') && (x(1) > 0);

    addRequired(inp, 'pos', checkPosDimensions);
    addRequired(inp, 'binWidth', checkBinWidth);
    addRequired(inp, 'shape', checkShape);
    addRequired(inp, 'dimensions', checkDims);

    parse(inp, pos, binWidth, shape, dimensions);

    t = pos(:, 1);
    x = pos(:, 2);
    if size(pos, 2) > 2
        y = pos(:, 3);
    else
        y = [];
    end

    binWidthX = binWidth(1);
    if length(binWidth) == 1
        binWidthY = binWidthX;
    else
        binWidthY = binWidth(2);
    end

    % filter out points that lie outside arena. This is important for circles
    % because they do not fully cover all the bins.
    if shape == bntConstants.ArenaShape.Circle
        radius = dimensions(1)/2;
        ind = sqrt(x.^2 + y.^2) > radius;
        x(ind) = nan;
        y(ind) = nan;
    end

    limitsX = [-dimensions(1)/2 dimensions(1)/2];
    [xBinned, nBinsX, edgesX] = helpers.bin(x, limitsX, binWidthX);
    nBins = nBinsX;

    if ~isempty(y)
        if length(dimensions) == 1
            if shape == bntConstants.ArenaShape.Track
                yLength = 1;
            else
                yLength = dimensions(1);
            end
        else
            if shape == bntConstants.ArenaShape.Track && dimensions(2) == 0
                yLength = 1;
            else
                yLength = dimensions(2);
            end
        end

        limitsY = [-yLength/2 yLength/2];
        [yBinned, nBinsY, edgesY] = helpers.bin(y, limitsY, binWidthY);
        nBins = [nBinsX nBinsY];
    end

    dt = diff(t);
    dt(end+1) = dt(end);

    if isempty(y)
        occupancy = general.accumulate(xBinned, dt, nBinsX)';
    else
        occupancy = general.accumulate([xBinned yBinned], dt, nBins)';
    end

    switch shape
        case bntConstants.ArenaShape.Circle
            radius = dimensions(1)/2;
            radiusBin = ceil(radius / binWidthX) + 1;

            halfSize = ceil(size(occupancy)/2);
            [rr, cc] = meshgrid(1:nBinsX, 1:nBinsY);

            distMap = sqrt((cc - halfSize(2)).^2 + (rr - halfSize(1)).^2); % each element is the distance
                                                % from the middle of the map to current point
            outerCircle = distMap >= radiusBin;
            distMap(outerCircle) = nan;

            numBins = sum(sum(isfinite(distMap)));
            coverage = (length(find(occupancy > 0)) / numBins) * 100;
        otherwise
            coverage = length(find(occupancy > 0)) / prod(nBins) * 100;
    end
end
