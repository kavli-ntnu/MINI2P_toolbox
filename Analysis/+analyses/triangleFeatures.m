% Calculate feature vector for a triangle
%

% Feature vector is:
%   centroid_x
%   centroid_y
%   orientation_1
%   orientation_2
%   orientation_3
%   area
function features = triangleFeatures(points)
    features = zeros(6, 1);

    fCentrX = 1;
    fCentrY = 2;
    fOr1 = 3;
    fOr2 = 4;
    fOr3 = 5;
    fArea = 6;

    centr = centroid(points);
    features(fCentrX) = centr(1);
    features(fCentrY) = centr(2);

    line = createLine(centr, points(1, :));
    features(fOr1) = rad2deg(lineAngle(line));

    line = createLine(centr, points(2, :));
    features(fOr2) = rad2deg(lineAngle(line));

    line = createLine(centr, points(3, :));
    features(fOr3) = rad2deg(lineAngle(line));
    
%     features(fArea) = triangleArea(points(1, :), points(2, :), points(3, :));

    features(6) = points(1, 1);
    features(7) = points(1, 2);
    features(8) = points(2, 1);
    features(9) = points(2, 2);
    features(10) = points(3, 1);
    features(11) = points(3, 2);
end
