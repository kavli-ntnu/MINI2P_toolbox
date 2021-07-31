% Calculate maximum or minimum of given function using quadratic interpolation
%
% This function first interpolates the given function around it's maximum/minimum
% value, then the maximum/minimum is calculated as vertex of quadratic function.
% This function is intended to be used for obtaining a maximum/minimum of a
% binned function. After binning, the maximum/minimum value will be given with
% bin resolution. By using the interpolation, it is possible to increase precision.
% Units of vx/vy are determined by units of bins/y.
%
%  USAGE
%   [vx, vy] = quadraticInterpolation(bins, y, limits, <useMax>)
%   bins        Bin/x values for function y.
%   y           Function for which maximum value is obtained.
%   limits      [l1 l2] two element vector that defines area of interpolation.
%               The function y is interpolated on interval [max(y)-l1 min(y)+l2].
%   <useMax>    Optional, boolean flag. If useMax is true, interpolation uses
%               max value of y. When useMax is false, interpolation uses min
%               value of y.
%   vx          X-value of maximum of y. It can be NaN if quadratic function has
%               it's minimum at that point instead of maximum, or if vx goes beyond
%               the range of bins.
%   vy          Y-value of maximum of y. It can be NaN if vx is NaN.
%
function [vx, vy] = quadraticInterpolation(bins, y, limits, useMax)
    if isempty(y)
        vx = nan;
        vy = nan;
        return;
    end
    if nargin < 4
        useMax = false;
    end
    limits = abs(limits);
    if useMax
        [~, i] = max(y);
    else
        [~, i] = min(y);
    end
    self = max(1, i-limits(1)):min(length(bins), i+limits(2));
    p = polyfit(bins(self), y(self), 2);
    vx = -p(2)/(2*p(1)); % x-value of quadratic vertex, standard equation

    % assign NaN to data with maxima beyond the bin range. In some instances
    % function defined by p defines a parabola with the contrary sign expected
    % and a vertex within bins (for POS speed cells p(1) should be negative and
    % for for NEG speed cells p(1) should be positive). Those values aren't
    % correct and are replaced by NaN.
    cond(1) = vx > bins(end) || vx < bins(1);
    cond(2) = useMax && p(1) > 0;
    cond(3) = ~useMax && p(1) < 0;
    if any(cond)
        vx = nan;
    end
    vy = polyval(p, vx);
end
