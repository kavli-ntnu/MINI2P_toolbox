% Assumes depth == 0
function point3D = proj2Dto3D(point2D, cameraParams, rotationMatrix, translationVector)
    N = size(point2D, 1);
    if N ~= 2
        point2D = point2D';
    end
    if size(point2D, 1) == 2
        point2D(3, :) = 1;
    end
    numPoints = size(point2D, 2);

    % if length(point2D) == 2
    %     point2D(3) = 1;
    % end
    % if size(point2D, 2) > 1
    %     point2D = point2D';
    % end
    K = cameraParams.IntrinsicMatrix';
    z_small = inv(rotationMatrix) * inv(K) * point2D;
    lamda = repmat(-translationVector(3), 1, numPoints) ./ z_small(3, :);

    translationVector = translationVector';
    % multiply lamda*z_small, i.e. lamda(i)*z_small(:, i)
    B = cellfun(@(x, y)(translationVector + x * y), num2cell(lamda), num2cell(z_small, 1), 'uniformoutput', false);
    point3D = cell2mat(B);
end