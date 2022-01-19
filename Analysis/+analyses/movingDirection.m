% Calculate moving direction of an animal
%
% Calculates moving direction for each position sample. Direction
% is defined as east = 0 degrees, north = 90 degrees, west = 180 degrees,
% south = 270 degrees. Direction is set to NaN for missing samples.
% Moving direction can be used to show that animal had no preferences in
% direction during running. Or as an analogue for occupancy map.
%
%  USAGE
%   [md, newPos] = analyses.movingDirection(pos, <options>)
%   pos         Matrix of position samples. Could be of size Nx3 or Nx5.
%   <options>   optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties         Values
%   -------------------------------------------------------------------------
%    'windowPoints'     [n1 n2], integers. The moving direction for current sample
%                       is calculated using the mean value of it's neighbours.
%                       windowPoints specifies how many points will be taken before
%                       the current sample (n1) and how many after (n2). Default
%                       value is [1 1], which means that the moving direction for
%                       sample i is determined by samples i-1 and i+1.
%                       Note that [1 1] gives you a histogram with very sharp
%                       bins at [0 90 180 270] degrees.
%    'step'             Integer. Moving direction can be calculated not for every
%                       position sample, but for every 'step' sample.
%                       Default is 1.
%   =========================================================================
%
%   md          Moving directions in degrees. If size of pos is Nx3, then
%               size of md is Nx1. If size of pos is Nx5, then size of md is Nx2.
%   newPos      Specifying 'step' greater than 1 has the same effect on position
%               data as changing the sampling rate. For example length(md) with
%               'step' == 2 should be N/2. But keeping the original number of position
%               samples is important because we match spikes to positions.
%               newPos contains a copy of input pos matrix with samples that are
%               excluded by moving direction calculation set to NaN. If 'step' == 2,
%               this means that every second sample in newPos will be set to NaN.
%
function [md, newPos] = movingDirection(pos, varargin)
    newPos = pos;
    inp = inputParser;
    defaultPoints = [1 1];
    defaultStep = 1;

    checkPosDimensions = @(x) size(x, 2) == 3 || size(x, 2) == 5;
    checkWindowPoints = @(x) (length(x) == 1 || length(x) == 2) && helpers.isivector(x, '>=1');
    checkStep = @(x) length(x) == 1 && helpers.isivector(x, '>=1');


    addRequired(inp, 'pos', checkPosDimensions);
    addParamValue(inp, 'windowPoints', defaultPoints, checkWindowPoints);
    addParamValue(inp, 'step', defaultStep, checkStep);

    parse(inp, pos, varargin{:});

    nBefore = inp.Results.windowPoints(1);
    if length(inp.Results.windowPoints) == 1
        nAfter = nBefore;
    else
        nAfter = inp.Results.windowPoints(2);
    end
    step = inp.Results.step;

    if isnan(nBefore) || isnan(nAfter) || isnan(step)
        error('MATLAB:InputParser:ArgumentFailedValidation', 'Either ''windowPoints'' or ''step'' contains NaN values. This is not supported');
    end


    % Number of position samples
    numSamples = size(pos, 1);
    mdInd = nBefore+1:step:numSamples-nAfter;

    kernel = [ones(1, nBefore)/nBefore 0 -ones(1, nAfter)/nAfter]; % minus in second part is not important

    droppedSamples = setdiff(1:size(pos, 1), mdInd);
    newPos(droppedSamples, 2:end) = nan;

    leftInd = mdInd - 1;
    rightInd = mdInd + 1;

    if size(pos, 2) > 3
        md = nan(numSamples, 2);
        md(mdInd, 1) = calcDirection(pos(:, bntConstants.PosX), pos(:, bntConstants.PosY), ...
            kernel, mdInd);
        md(mdInd, 2) = calcDirection(pos(:, bntConstants.PosX2), pos(:, bntConstants.PosY2), ...
            kernel, mdInd);
    else
        md = nan(numSamples, 1);
        md(mdInd, 1) = calcDirection(pos(:, bntConstants.PosX), pos(:, bntConstants.PosY), ...
            kernel, mdInd);
    end
end

function md1 = calcDirection(x, y, kernel, ind)
    cr = conv(x, kernel, 'same');
    X = cr(ind);
    cr = conv(y, kernel, 'same');
    Y = cr(ind);
    md1 = mod(atan2d(Y, X), 360);
end
