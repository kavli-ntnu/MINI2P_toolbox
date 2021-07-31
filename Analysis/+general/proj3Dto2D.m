% Project a 3D point to 2D using a pinhole camera model
%
% Projection is done according to Eq. 2.18 from
%
%  USAGE
%   point2D = general.proj3Dto2D(point3D, cameraParams, rotationMatrix, translationVector)
%   point3D             Nx3 set of coordinates for points. Column-based matrix in form [x y z].
%   cameraParams        cameraParameters object. Object for storing camera parameters,
%                       specified as a cameraParameters returned by the estimateCameraParameters
%                       function or the Camera Calibrator app. This object contains the intrinsic,
%                       extrinsic, and lens distortion parameters of a camera.
%   rotationMatrix      3x3 matrix describing camera rotation. Output of extrinsics.
%   translationVector   1x3 vector describing camera translation. Output of extrinsics.
%   point2D             Nx2 projected points. Matrix is organized column-based in form [x y].
%
function point2D = proj3Dto2D(point3D, cameraParams, rotationMatrix, translationVector)
    K = cameraParams.IntrinsicMatrix';

    N = size(point3D, 1);
    if N ~= 3
        point3D = point3D';
    end
    numPoints = size(point3D, 2);
    
    point2D = K*rotationMatrix*point3D - repmat((K*rotationMatrix*translationVector'), 1, numPoints);
    
    % Divide each column by the last element in that column
    B = cellfun(@(x)(x ./ x(end)), num2cell(point2D, 1), 'uniformoutput', false);
    point2D = cell2mat(B);
    
    point2D(3, :) = [];
    point2D = point2D';
end