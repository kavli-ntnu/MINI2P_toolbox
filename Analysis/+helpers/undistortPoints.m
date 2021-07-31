% Undistort points using camera parameters
%
% This function is much faster than Matlab's undistortPoints and works on large
% point sets. The code is based on cvUndistortPoints function from OpenCV library.
% Also this function works with NaN values.
%
%  USAGE
%   undistortedPoints = undistortPoints(points, cameraParams)
%   points              Nx2 ([x y]) matrix of points.
%   cameraParams        Camera calibration information.
%   undistortedPoints   Nx2 ([x y]) matrix with undistorted points.
%
function undistortedPoints = undistortPoints(points, cameraParams)
    origin = [cameraParams.IntrinsicMatrix(3, 1:2)];
    cx = origin(1);
    cy = origin(2);
    fx = cameraParams.IntrinsicMatrix(1, 1);
    fy = cameraParams.IntrinsicMatrix(2, 2);
    ifx = 1 / fx;
    ify = 1 / fy;
    iters = 50;
    k = cameraParams.RadialDistortion;
    p_coef = cameraParams.TangentialDistortion;
    
    if length(k) == 2
        k(3) = 0;
    end

    x = points(:, 1);
    y = points(:, 2);
    x = (x - cx) .* ifx;
    y = (y - cy) .* ify;
    x0 = x;
    y0 = y;

    for j = 1:iters
        r2 = x.*x + y.*y;
        icdist = 1 ./ (1 + ((k(3).*r2 + k(2)).*r2 + k(1)) .* r2);
        deltaX = 2.*p_coef(1).*x.*y + p_coef(2).*(r2 + 2.*x.*x);
        deltaY = p_coef(1) .* (r2 + 2.*y.*y) + 2.*p_coef(2).*x.*y;
        x = (x0 - deltaX) .* icdist;
        y = (y0 - deltaY) .* icdist;
    end

    undistortedPoints(:, 1) = round(x .* fx + cx);
    undistortedPoints(:, 2) = round(y .* fy + cy);
end