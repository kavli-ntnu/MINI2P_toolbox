% Smooth 1D data with a moving average Gaussian kernel
%
% The function creates a Gaussian kernel with specified standard deviations
% and applies it to the underlying data.
% It differs from general.smooth function in that it is possible to specify
% desired width of the Gaussian. Also, this new function has been introduced
% for backward compatibility, i.e. smoothing of rate maps and turning curve
% is still done by general.smooth.
%
%  USAGE
%   smoothed = general.smoothGauss(data, sigma, <len>)
%   data    data to smooth, 1D vector.
%   sigma   standard deviation for Gaussian kernel, measured
%           in number of samples (0 = no smoothing).
%   len     Optional length of Gaussian filter. If not specified,
%           then 6*sigma+1 is used.
%
function smoothed = smoothGauss(data, sigma, len)
    if ~isvector(data)
        error('Smoothing applies only to vectors (type ''help <a href="matlab:help general.smoothGauss">general.smoothGauss</a>'' for details).');
    end
    % Vectors must be 'vertical'
    if size(data, 1) == 1
        data = data';
    end
    if sigma <= 0
        smoothed = data;
        return;
    end
    if nargin < 3
        len = 6*sigma + 1;
    end

    % define gaussian filter
    paramAlpha = (len - 1) / (2 * sigma);
    gw = gausswin(len, paramAlpha);
    gw = gw / sum(gw);

    % apply filter
    smoothed = conv(data, gw, 'same');
end