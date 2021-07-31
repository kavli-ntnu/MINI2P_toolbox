% Generate 2D image with Gaussian
% mu                [x y] pairs.
% sigma             Nx2 matrix of sigmas for x and y dimension.
% outputSize        Size of the output matrix, one or two elements. This is direct
%                   input to zeros() function.
% amplitude         Optional vector of amplitudes. Default is all 1.
%
function res = gauss2d(mu, sigma, outputSize, amplitude)
    numPoints = size(mu, 1);
    numSigmas = size(sigma, 1);
    if numPoints ~= numSigmas
        error('Invalid size');
    end
    if nargin < 4
        amplitude = ones(numPoints, 1);
    end

    res = zeros(outputSize);
    gsize = size(res);
    [R, C] = ndgrid(1:gsize(1), 1:gsize(2));
    for i = 1:numPoints
        if isempty(find(isnan(mu(i, :)), 1))
            res = res + gaussC(C, R, sigma(i, :), mu(i, :), amplitude(i));
%             plot.colorMap(res);
        end
    end
end

function val = gaussC(x, y, sigma, center, amplitude)
    xc = center(1);
    yc = center(2);
    sigma_x = sigma(1);
    if length(sigma) > 1
        sigma_y = sigma(2);
    else
        sigma_y = sigma_x;
    end

    % exponent = ((x - xc).^2 + (y-yc).^2) ./ (2*sigma(1));

    exponent = ( (((x - xc).^2)./(2*sigma_x)) + (((y-yc).^2) ./ (2*sigma_y)));

    val = amplitude * exp(-exponent);
end