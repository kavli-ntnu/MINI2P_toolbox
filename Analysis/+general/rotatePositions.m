% Rotate position samples
%
% This function rotates provided position samples to specified angle.
%
%  USAGE
%   rpos = rotatePositions(pos, degAngle)
%   pos         Matrix with position samples.
%   degAngle    Rotation angle in degrees.
%   rpos        Rotated position samples. Same dimension as pos.
%
function rpos = rotatePositions(pos, degAngle)
    rpos(:, 2) = pos(:, 2)*cosd(degAngle) - pos(:, 3)*sind(degAngle);
    rpos(:, 3) = pos(:, 2)*sind(degAngle) + pos(:, 3)*cosd(degAngle);
    if size(pos, 2) > 3
        rpos(:, 4) = pos(:, 4)*cosd(degAngle) - pos(:, 5)*sind(degAngle);
        rpos(:, 5) = pos(:, 4)*sind(degAngle) + pos(:, 5)*cosd(degAngle);
    end
end