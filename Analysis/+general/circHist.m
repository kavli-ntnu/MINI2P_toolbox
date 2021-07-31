% Historgram count for circular data
%
% Creates histogram of circular data. It is a wrapper around <a href="matlab:help histc">histc</a> function.
%
%  USAGE
%   [n, edges, midPoints] = general.circHist(theta, binWidth)
%   theta           Vector of angles in degrees.
%   binWidth        Width of each bin expressed in degrees. Bins are spaced in the range [0, 360].
%   n               Vector that contains the histogram count values. Number of elements is determined
%                   based on binWidth parameter.
%   edges           Bins edge values, length(edges) = lenth(n) + 1. Values theta(i) are assigned to bin k
%                   if edges(k) <= theta(i) < edges(k+1).
%   midPoints       Vector of middle point values for each bin. If edge(1) == 30 and edge(2) == 60, then
%                   midPoints(1) == 45. This value is usefull for further processing of vector n.
%                   For example one should use it to calculate mean angle. Units are degrees.
%
%  SEE
%   See also histc, circ_mean
%
function [n, edges, midPoints] = circHist(theta, binWidth)
    twoPi = 360;
    theta = mod(theta, twoPi); % Make sure 0 <= theta <= 2*pi

    numBins = ceil(twoPi / binWidth);
    x = (0:numBins-1) * twoPi/numBins + 180./numBins;
    midPoints = x;
%     edges = sort(mod([(x(2:end) + x(1:end-1))/2 (x(end) + x(1) + twoPi)/2], twoPi));
%     edges = [edges edges(1) + twoPi];

    [n, edges] = histcounts(theta, numBins, 'binlimits', [0 twoPi]);
%     if ~isempty(n)
%         n(end-1) = n(end-1) + n(end);
%         n(end) = [];
%     end
end