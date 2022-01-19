% Detect left and right runs on a linear track
%
% This function locates left and right runs of an animal on a linear track.
% The runs are located by finding peak values of animal's x-coordinate data. Peaks determines
% situations when an animal have been running in one direction and then turned around and run
% in the opposite direction.
%
%  USAGE
%   runData = analyses.findRuns(pos, <options>)
%   pos            Nx3 or Nx5 matrix, position data
%   <options>      optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties    Values
%   -------------------------------------------------------------------------
%    'borderThreshold'      Defines the search zone around the border of a track.
%                           The value should be given in range 0.5-1. It corresponds to
%                           the percentage of the track length. For example, if track length is 100 cm
%                           and borderThreshold is 0.6 (track is zero centred), then the function
%                           will look for peaks in two zones: 30..50 cm and -30..-50 cm. Useful
%                           if animal doesn't run until the very end of the track.
%                           Default value is 0.6.
%    'omitEdgeActivity'     Sometimes animal doesn't immediately run back on the track, rather
%                           stays at the end/beginning of the track for some time. If you don't
%                           want to include such activities in the final left/right runs, then set
%                           this property to 'on'. Otherwise set it to 'off'. Default value is 'on'.
%    'minDuration'          Value of a minimum run duration in seconds. If provided, then runs that are
%                           shorter than this value, are discarded. Default value is 0, meaning all
%                           runs are kept.
%    'minLength'            Value of a minimum run length in cm. If provided, then runs that
%                           are shorter than this value will be discarded. Default value is 0, meaning
%                           that all runs are kept.
%    'selectivity'          Selectivity of peak finder algorithm. See parameter
%                           sel of function peakfinder. Default is 1!
%   -------------------------------------------------------------------------
%   runData         Array of structures with fields:
%               indices     [runStartPos runEndPos]
%               isLeftRun   TRUE if it is left run, FALSE if it is right run
%
function runData = findRuns(pos, varargin)
    % a small hack to handle Virmen recordings.
    global gBntData;
    global gCurrentTrial;

    inp = inputParser;

    defaultBorderThreshold = 0.6;
    defaultOmitEdgeActivity = 'on';
    defaultMinDuration = 0;
    defaultSelectivity = 1;
    defaultMinLength = 0;

    checkPosDimensions = @(x) size(x, 2) >= 3;
    checkBorderThreshold = @(x) helpers.isdscalar(x, '>=0.5', '<=1');
    checkOmitEdge = @(x) strcmpi(x, 'on') || strcmpi(x, 'off');
    checkSelectivity = @(x) isempty(x) || helpers.isdscalar(x, '>0');
    checkDuration = @(x) helpers.isdscalar(x, '>=0');

    addRequired(inp, 'pos', checkPosDimensions);
    addParamValue(inp, 'borderThreshold', defaultBorderThreshold, checkBorderThreshold);
    addParamValue(inp, 'omitEdgeActivity', defaultOmitEdgeActivity, checkOmitEdge);
    addParamValue(inp, 'minDuration', defaultMinDuration, checkDuration);
    addParamValue(inp, 'selectivity', defaultSelectivity, checkSelectivity);
    addParamValue(inp, 'minLength', defaultMinLength, checkDuration);

    parse(inp, pos, varargin{:});

    % get parsed results
    borderThreshold = inp.Results.borderThreshold;
    omitEdgeActivity = strcmpi(inp.Results.omitEdgeActivity, 'on');
    duration = inp.Results.minDuration;
    selectivity = inp.Results.selectivity;
    minLength = inp.Results.minLength;

    posx = pos(:, 2);

    if ~helpers.isstring(gBntData{gCurrentTrial}.system, bntConstants.RecSystem.Virmen)
        if length(posx) < 100000
            % otherwise it will never end
            span = round(length(posx) * 0.015);
            y = smooth(posx, span, 'loess');
        else
            y = medfilt1(posx, 15);
            y = medfilt1(y, 15);
        end
    else
        y = medfilt1(posx, 15);
        y = medfilt1(y, 15);
    end
    upperThreshold = nanmax(y) * borderThreshold;
    lowerThreshold = nanmin(y) * borderThreshold;

    peaks = peakfinder(y, selectivity); % 1 because we want to detect peaks around track ends
    indToDel = false(1, length(peaks));
    for i = 1:length(peaks)
        curPeak = peaks(i);
        if y(curPeak) < upperThreshold && posx(curPeak) < upperThreshold
            indToDel(i) = true;
            continue;
        end
    end
    peaks(indToDel) = [];

    peaksn = peakfinder(-y, selectivity);
    indToDel = false(1, length(peaksn));
    for i = 1:length(peaksn)
        curPeak = peaksn(i);
        if y(curPeak) > lowerThreshold && posx(curPeak) > lowerThreshold
            indToDel(i) = true;
            continue;
        end
    end
    peaksn(indToDel) = [];

    peaks = sort(cat(1, peaks, peaksn));

    % let's calculate signal's Fs by looking somewhere in the middle
    ind = round(length(posx) / 2);
    Fs = 1 / (pos(ind, 1) - pos(ind-1, 1));

    runData = [];
    for i = 1:length(peaks)
        if i == 1
            beginPoint = 1;
        else
            beginPoint = peaks(i-1);
        end
        s = std(y(beginPoint:peaks(i)));
        if s < 5 && omitEdgeActivity % 5 is just an arbitary threshold. Perhaps it should depend
                                     % on track length. Run on a 100 cm track produces STD of ~25-30.
            continue;
        end
        tLength = peaks(i) - beginPoint;
        length_sec = tLength / Fs;
        if length_sec < duration
            continue;
        end
        length_cm = abs(posx(peaks(i)) - posx(beginPoint));
        if length_cm < minLength
            continue;
        end

        runData = [runData; struct('indices', [beginPoint peaks(i)], 'isLeftRun', false)];
        if y(beginPoint) > y(peaks(i))
            runData(end).isLeftRun = false;
        else
            runData(end).isLeftRun = true;
        end
    end
end

