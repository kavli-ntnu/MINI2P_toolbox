% Calculate statistics of a turning curve.
%
% Calculates various statistics for a turning curve.
% 1. Mean vector length of a head direction rate map.
% The value will range from 0 to 1. 0 means that there are so much dispersion
% that a mean angle cannot be described. 1 means that all data are
% concentrated at the same direction. Note that 0 does not necessarily
% indicate a uniform distribution.
% Calculation is based on Section 26.4, J.H Zar - Biostatistical Analysis 5th edition,
% see eq. 26.13, 26.14.
%
%  USAGE
%
%    tcStat = analyses.tcStatistics(turningCurve, binWidth, percentile)
%    turningCurve       Data of a turning curve. It can be a vector or a matrix.
%                       In case of matrix only the second column is used.
%    binWidth           Bin width of turning curves in degrees.
%    percentile         Percentile value for the head direction arc calculation.
%                       Arc is between two points with values around
%                       globalPeak * (percentile / 100). Value should be in
%                       range 0..100.
%    tcStat             Structure with various statistics:
%       mean            Mean direction in degrees in range 0..360.
%       r               Mean vector length.
%       std             Circular standard deviation in radians.
%       score           Head direction score.
%       peakRate        Turning curve peak rate.
%       meanRate        Turning curve mean rate.
%       peakDirection   Turning curve peak direction in degrees.
%       arcAngRad       Angle in radians of arc around the peak. Arc is between two points that lie on
%                       globalPeakValue * (percentile / 100) distance away of global peak.
%       halfCwInd       Bin number that corresponds to the 'right' edge of the arc around the peak.
%       halfCwIndRad    Value in radians that correspond to the angle of the right edge of the arc.
%       halfCcwInd      Bin number of the 'left' edge of the arc around the peak.
%       halfCcwIndRad   Value in radians that correspond to the angle of the left edge of the arc.
%
%
%  NOTE
%   1. This function assumes that the data is grouped (which is true for a turning curve).
%      You should not use this function if you want to calculate statistics of a set of angles.
%
function tcStat = tcStatistics(turningCurve, binWidth, percentile)

    tcStat.mean = nan;
    tcStat.r = nan;
    tcStat.std = nan;
    tcStat.score = nan;
    tcStat.peakRate = nan;
    tcStat.peakDirection = nan;
    tcStat.meanRate = nan;
    tcStat.arcAngRad = nan;
    tcStat.halfCwInd = nan;
    tcStat.halfCwIndRad = nan;
    tcStat.halfCcwInd = nan;
    tcStat.halfCcwIndRad = nan;

    if isempty(turningCurve)
        return;
    end

    if size(turningCurve, 2) > 1
        % we have mapAxis
        mapAxis = CircStat2012a.circ_ang2rad(turningCurve(:, 1));
        turningCurve = turningCurve(:, 2);
    else
        numBins = length(turningCurve);
        binWidth = CircStat2012a.circ_ang2rad(ceil(360 / numBins));
        halfBin = binWidth / 2;

        mapAxis = halfBin:binWidth:(numBins*binWidth)-halfBin;
        mapAxis = mapAxis';
    end

    if size(turningCurve, 2) > 2
        % transpose to make the same dimensions as axis
        turningCurve = turningCurve';
    end
    
    % remove nans
    nans = isnan(turningCurve);
    turningCurve(nans) = [];
    mapAxis(nans) = [];
    
    if isempty(turningCurve)
        return;
    end

    tcStat.mean = mod(CircStat2012a.circ_rad2ang(CircStat2012a.circ_mean(mapAxis, turningCurve)), 360);
    tcStat.r = CircStat2012a.circ_r(mapAxis, turningCurve);
    tcStat.std = sqrt(2 * (1 - tcStat.r)); % Eq. 26.20 from J. H. Zar

    [globalPeak, globalPeakInd] = nanmax(turningCurve);
    halfPeak = globalPeak * (percentile / 100);
    tcStat.halfCwInd = find(turningCurve(1:globalPeakInd) <= (globalPeak - halfPeak), 1, 'last');
    if isempty(tcStat.halfCwInd)
        % try another search for case when peak is around 0 deg.
        tcStat.halfCwInd = find(turningCurve(globalPeakInd+1:end) <= (globalPeak - halfPeak), 1, 'last');
        if isempty(tcStat.halfCwInd)
            tcStat.halfCwInd = 0;
        else
            % we have counter-clockwise direction
            tcStat.halfCwInd = globalPeakInd + tcStat.halfCwInd;
        end
        tcStat.halfCwIndRad = deg2rad(tcStat.halfCwInd * binWidth);
    else
        tcStat.halfCwIndRad = deg2rad(tcStat.halfCwInd * binWidth);
    end

    tcStat.halfCcwInd = find(turningCurve(globalPeakInd+1:end) <= (globalPeak - halfPeak), 1, 'first');
    if isempty(tcStat.halfCcwInd)
        tcStat.halfCcwInd = find(turningCurve(1:globalPeakInd) <= (globalPeak - halfPeak), 1, 'first');
        if isempty(tcStat.halfCcwInd)
            tcStat.halfCcwInd = length(turningCurve);
            tcStat.halfCcwIndRad = deg2rad(tcStat.halfCcwInd * binWidth);
        else
            tcStat.halfCcwIndRad = deg2rad(tcStat.halfCcwInd * binWidth) + (2*pi);
        end
    else
        tcStat.halfCcwInd = globalPeakInd + tcStat.halfCcwInd;
        tcStat.halfCcwIndRad = deg2rad(tcStat.halfCcwInd * binWidth);
    end

    % arc angle in radians
    tcStat.arcAngleRad = tcStat.halfCcwIndRad - tcStat.halfCwIndRad;

    tcStat.score = 1 - tcStat.arcAngleRad / pi;

    tcStat.peakRate = globalPeak;
    tcStat.meanRate = nanmean(turningCurve);
    tcStat.peakDirection = CircStat2012a.circ_rad2ang(mapAxis(globalPeakInd));
end
