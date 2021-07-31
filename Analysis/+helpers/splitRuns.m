% Split position data into left and right runs
%
function [leftPos, rightPos] = splitRuns(pos, leftIndices, rightIndices)

    leftPos = nan(size(pos));
    leftPos(:, 1) = pos(:, 1);
    for j = 1:size(leftIndices, 1)
        leftPos(leftIndices(j, 1):leftIndices(j, 2), :) = pos(leftIndices(j, 1):leftIndices(j, 2), :);
    end

    rightPos = nan(size(pos));
    rightPos(:, 1) = pos(:, 1);
    for j = 1:size(rightIndices, 1)
        rightPos(rightIndices(j, 1):rightIndices(j, 2), :) = pos(rightIndices(j, 1):rightIndices(j, 2), :);
    end
end