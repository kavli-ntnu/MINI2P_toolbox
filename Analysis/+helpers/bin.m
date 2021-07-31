% Wrapper around Matlab's histc function with some modifications
%
% This function assigns each input value a bin number according to specified
% limits of values and width of each bin.
% Vector of bin edges is created out of specified limits and bin width. Value
% x(i) is assigned to bin k if edges(k) <= x(i) < edges(k+1). So far this is
% exactly what Matlab's function histc does. The difference comes in the last bin.
% Condition of assignment for the last bin is edges(end-1) <= x(i) <= edges(end).
%
%  USAGE
%   [bins, nBins, edges] = helpers.bin(x, limits, binWidth)
%   x           1D vector of values that should be binned.
%   limits      [min_x max_x] limit values that define bins. These need not be equal to
%               min(x) and max(x). Values of x that equal to min_x or max_x are all
%               included in the final histogram.
%   binWidth    width of each individual bin. Used to create edges vector.
%   bins        Vector of assigned bin numbers to x. length(bins) == length(x).
%               Values outside of limits are assigned to NaN.
%   nBins       Number of bins.
%   edges       Edge values of bins.
%
function [bins, nBins, edges] = bin(x, limits, binWidth)
    nBins = ceil((limits(2) - limits(1)) / binWidth);

    edges = limits(1):binWidth:limits(2);
    if length(edges) == nBins
        edges(end+1) = limits(2);
    end

    ind = (x == edges(end)); % find indices of values that according to mathematics should be
                                 % assigned to bin nBins+1
    [~, bins] = histc(x, edges);
    bins(ind) = nBins;
    bins(bins > nBins) = nan;
    bins(bins == 0) = nan;
end
