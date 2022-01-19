% Detect left and right runs on a linear track using threshold method
%
% This function detects runs based on position threshold. Two position points
% are calculated for the beginning and end of the track (called thresholds). Whenever
% animal reaches threshold, information about the track is stored.
%
%  NOTE
%   The function assumes that position data is zero centered! This means that the
%   range of x-values shall be within [-trackLength/2 trackLength/2].
%
%  USAGE
%    [leftIndices, rightIndices] = analyses.findRunsByThreshold(pos, thresholdLevels, <options>)
%    pos                Nx3 or Nx5 matrix, position data. Position data can be either
%                       zero centerd (x values in range [-trackLength/2 trackLength/2])
%                       or not (x values in range [0 trackLength]).
%    thresholdLevels    [v1] or [v1 v2]. Threshold values that will be applied to pos directly.
%                       This means that units of thresholdLevels should be the same as pos units.
%                       See also general.runThreshold.
%    <options>          Optional list of property-value pairs (see table below)
%
%    ==============================================================================================
%     Properties        Values
%    ----------------------------------------------------------------------------------------------
%     'minSamples'      Minimum number of samples in a run. If run contains less samples, then it
%                       is discarded. In other words, this argument controls minimum duration of a run.
%                       Default value is 0, meaning that all runs are kept.
%     'minLength'       Minimum length of a run given in the same units as pos (i.e. in cm). Runs
%                       of shorter length than minLength will be discarded.
%                       Default value is 0, meaning that all runs are kept.
%   ==============================================================================================
%
%    leftIndices        Nx1 vector that indicates runs in pos matrix. The vector will have 1 values for
%                       position samples that correspond to the first run, 2 values for position
%                       samples that corrrespoind to the second run and so on. This can be usefull
%                       to quickly filter the data. For example, assign NaN values for everything
%                       that is not a run:
%                       pos(leftIndices == 0, :) = nan;
%                       If no left runs have been detected, then max(leftIndices) is 0.
%    rightIndices       Same as leftIndices, but for the right runs.

%    leftIndices        Lx2 matrix. L - number of detected left runs.
%                       Columns contain row numbers of pos matrix, which mark
%                       the start and the end of the run. First column is for start.
%                       Second column is for the end. For example, matrix
%                       1   123
%                       567 1115
%                       describes to left runs. First run starts at position with index 1
%                       and stops at position with index 123. Second run starts at 567
%                       and ends at 1115.
%                       The matrix is empty if no left runs has been detected.
%    rightIndices       Rx2 matrix. R - number of detected right runs. It has the same structure
%                       as leftIndices matrix.
%
function [leftIndices, rightIndices] = findRunsByThreshold(pos, thresholdLevels, varargin)
    inp = inputParser;

    defaultMinSamples = 0;
    defaultMinLength = 0;
    checkPosDimensions = @(x) size(x, 2) >= 2;
    checkThreshold = @(x) ~isempty(x) && length(x) == 2;

    addRequired(inp, 'pos', checkPosDimensions);
    addRequired(inp, 'thresholdLevels', checkThreshold);
    addParameter(inp, 'minSamples', defaultMinSamples, @(x) x >= 0);
    addParameter(inp, 'minLength', defaultMinLength, @(x) x >= 0);

    parse(inp, pos, thresholdLevels, varargin{:});

    lowerThresh = thresholdLevels(1);
    upperThresh = thresholdLevels(2);

    % prepare outputs
    leftIndices = zeros(size(pos, 1), 1);
    rightIndices = zeros(size(pos, 1), 1);

    posx = pos(:, bntConstants.PosX);

    isCurrentRunLeft = false;
    isRunning = false;
    curRunStart = 1;
    minSamples = inp.Results.minSamples;
    minLength = inp.Results.minLength;
    curLeftRun = 0; % counters
    curRightRun = 0;

    firstNonNan = find(~isnan(posx), 1);
    if posx(firstNonNan) < lowerThresh
        % animal is in lower illegal area.
        ind = find(posx > lowerThresh, 1); % find first position that lies inside the legal area
        isCurrentRunLeft = true;
        isRunning = true;
        curRunStart = ind;
    elseif posx(firstNonNan) >= upperThresh
        % animal is in upper illegal area
        ind = find(posx < upperThresh, 1);
        isCurrentRunLeft = false;
        isRunning = true;
        curRunStart = ind;
    end

    if ~isRunning
        % animal is not in the illegal zones, somewhere in between of the track.
        indLower = find(posx <= lowerThresh, 1);
        indUpper = find(posx >= upperThresh, 1);
        if indLower < indUpper
            % animal reached lower part first. In general now there should be
            % a run from lower part to the upper.
            curRunStart = indLower;
            isRunning = false; % animal is in illegal zone, not running yet
            nextRunIsLeft = true;
        else
            curRunStart = indUpper;
            isRunning = false;
            nextRunIsLeft = false;
        end
    end

    if isempty(curRunStart)
        return;
    end

    while true
        if isRunning
            if isCurrentRunLeft
                % -2 because we need -1 to get the actual correct index, another
                % -1 comes from the fact that we want the previous sample
                endInd = curRunStart + find(posx(curRunStart:end) >= upperThresh, 1) - 2;
                if isempty(endInd)
                    break;
                end

                isRunning = false;
                isCurrentRunLeft = false;
                nextRunIsLeft = false;

                potentialTurn = find(posx(curRunStart:endInd) <= lowerThresh, 1, 'last');
                if ~isempty(potentialTurn)
                    curRunStart = curRunStart + potentialTurn - 1 + 1;
                end

                if endInd - curRunStart >= minSamples && abs(posx(endInd) - posx(curRunStart)) > minLength
                    curLeftRun = curLeftRun + 1;
                    leftIndices(curRunStart:endInd) = curLeftRun;
                    %leftIndices = cat(1, leftIndices, [curRunStart endInd]);
                end

                curRunStart = endInd + 1;
            else
                % -2 because we need -1 to get the actual correct index, another
                % -1 comes from the fact that we want the previous sample
                endInd = curRunStart + find(posx(curRunStart:end) <= lowerThresh, 1) - 2;
                if isempty(endInd)
                    break;
                end

                potentialTurn = find(posx(curRunStart:endInd) >= upperThresh, 1, 'last');
                if ~isempty(potentialTurn)
                    curRunStart = curRunStart + potentialTurn - 1 + 1;
                end

                isRunning = false;
                nextRunIsLeft = true;

                if endInd - curRunStart >= minSamples && abs(posx(endInd) - posx(curRunStart)) > minLength
                    curRightRun = curRightRun + 1;
                    rightIndices(curRunStart:endInd) = curRightRun;
                    % rightIndices = cat(1, rightIndices, [curRunStart endInd]);
                end

                curRunStart = endInd + 1;
            end
            continue;
        end
        if nextRunIsLeft
            ind = find(posx(curRunStart:end) >= lowerThresh, 1);
            if isempty(ind)
                break;
            end
            curRunStart = ind + curRunStart - 1;
            isRunning = true;
            isCurrentRunLeft = true;
        else
            ind = find(posx(curRunStart:end) <= upperThresh, 1);
            if isempty(ind)
                break;
            end
            curRunStart = ind + curRunStart - 1;
            isRunning = true;
            isCurrentRunLeft = false;
        end
    end
end