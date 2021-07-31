% Calculate exact threshold levels for run detection on a linear track
%
% This function calculates the exact threshold levels (so that they can be
% directly applied to position data) based on the provided values.
% This function can be used to get the levels which can be further used
% to detect runs or/and to plot the thresholds.
% The threshold should be applied to both sides of the track, i.e. animal
% activity at the beginning of the track and at the end of the track should
% not be included in the analyses.
%
%  USAGE
%   th = general.runThreshold(trackLength, thresholdType, thresholdValue)
%   trackLength     Linear track length in whatever units your positions are.
%                   If you have positions in cm, then trackLength should also
%                   be in cm. No default value.
%   thresholdType   '%' or 'direct'. '%' means that the thresholdValue is given
%                   as the percentage of the trackLength. 'direct' means
%                   that thresholdValue has the same units as trackLength
%                   and it is the offset from the track edges.
%                   Default value is '%'
%   thresholdValue  [v] or [v1 v2]. Desired threshold value. If only a single value
%                   is given, then the same threshold level will be used for both
%                   sides of the track. If two values are provided, then v1 is used
%                   at the start of the track and v2 is used at the end of the track.
%                   The range of values depend on thresholdType. If thresholdType is '%',
%                   then thresholdValue should be given in range [0, 49].
%                   If thresholdType is 'direct', then thresholdValue should be within
%                   range [0, trackLength].
%                   Default value is 5 [%].
%   th              [lower upper] threshold levels. They are provided in range [0 trackLength].
%                   This means that if your positions are zero centred (there are negative values),
%                   then it is necessary to subtract trackLength/2 from th (th = th - trackLength/2).
%
function th = runThreshold(trackLength, thresholdType, thresholdValue)
    inp = inputParser;
    defaultType = '%';
    defaultValue = 5;

    checkType = @(x) strcmpi(x, '%') || strcmpi(x, 'direct');

    addRequired(inp, 'trackLength', @(x) ~isempty(x) && length(x) == 1 && x > 0);
    addOptional(inp, 'thresholdType', defaultType, checkType);
    addOptional(inp, 'thresholdValue', defaultValue, @(x) ~isempty(x) && length(x) <= 2);

    parse(inp, trackLength, thresholdType, thresholdValue);
    isPercent = strcmpi(thresholdType, '%');

    % additional check of thresholdValue
    if isPercent
        if any(thresholdValue < 0 | thresholdValue > 49)
            error('BNT:arg', 'Threshold type is percentage, and threshold value is out of range. Check description of thresholdValue argument.')
        end
    else
        if any(thresholdValue < 0 | thresholdValue > trackLength)
            error('BNT:arg', 'Threshold type is ''direct'', and threshold value is either negative or greater than the track length. Check your thresholdType value.');
        end
    end

    if isPercent
        thresholdValue = thresholdValue / 100;
        offset = thresholdValue * trackLength;
    else
        offset = thresholdValue;
    end

    if length(offset) == 1
        offset(2) = offset(1);
    end

    th(1) = offset(1);
    th(2) = trackLength - offset(2);
end