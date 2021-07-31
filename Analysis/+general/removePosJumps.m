% Remove bad tracking coordinates (position jumps)
%
% A position "jump" correspond to position samples that imply that the rat is
% moving quicker than physical possible.
%
%
function [x, y] = removePosJumps(x, y, threshold, stdThreshold)
    N = length(x);
    % Indexes to position samples that are to be removed
    remInd = zeros(N,1);
    remCounter = 0;

    diffX = diff(x);
    diffY = diff(y);
    diffR = sqrt(diffX.^2 + diffY.^2);
    ind = find(diffR > threshold);

    if isempty(ind)
        return;
    end

    if ind(end) == length(x)
        offset = 2;
    else
        offset = 1;
    end

    for ii = 1:length(ind)-offset
        if ind(ii+1) == ind(ii)+1
            % A single sample position jump, tracker jumps out one sample and
            % then jumps back to path on the next sample. Remove bad sample.
            remCounter = remCounter + 1;
            remInd(remCounter) = ind(ii)+1;
            ii = ii+1;
            continue
        else
            % Not a single jump. 2 possibilities:
            % 1. Tracker jumps out, and stay out at the same place for several
            % samples and then jumps back.
            % 2. Tracker just has a small jump before path continues as normal,
            % unknown reason for this. In latter case the samples are left
            % untouched.
            idx = find(x(ind(ii)+1:ind(ii+1)+1)==x(ind(ii)+1));
            if length(idx) == length(x(ind(ii)+1:ind(ii+1)+1));
                n = ind(ii+1)+1 - ind(ii);
                remInd(remCounter+1:remCounter+n) = (ind(ii)+1:ind(ii+1)+1)';
                remCounter = remCounter + n;
            end
        end
    end

    remInd = remInd(1:remCounter);

    % Remove the samples
    x(remInd) = NaN;
    y(remInd) = NaN;

    % there could be tracking outliers. They are commonly several values between longer list of
    % NaNs. These values lie far away from 'good' points, so discard everything that is
    % further than 2.5*STD
    med(1) = nanmin(x) + (nanmax(x) - nanmin(x))/2; % roughy middle point of the arena, better
    med(2) = nanmin(y) + (nanmax(y) - nanmin(y))/2; % than median or mean since it's less biased
    std = nanstd([x, y]);
    lowerBound = med - (stdThreshold*std);
    upperBound = med + (stdThreshold*std);

    x(x < lowerBound(1)) = nan;
    x(x > upperBound(1)) = nan;

    y(y < lowerBound(2)) = nan;
    y(y > upperBound(2)) = nan;
end