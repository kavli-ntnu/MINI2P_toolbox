% Fix point assignment to LEDs in tracked positions
%
% Axona doesn't use colour for tracking, but the size of LEDs. In case when
% only one LED is visible the assignment could be swapped. This means that
% points that should be treated as points of LED1 are recorded as points of
% LED2 (and vice versa). This function tries to detect this situations
% and reassign the points. It goes through all the points and checks how far
% are the current points from last tracked points of LED1 and LED2. Based on
% the distance, the function reassigns the current points.
% The function should be used in chain:
%   * it is better to fix isolated points prior to call of fixLedAssignment.
%   * it is also possible to interpolate positions in the beginning of the track if there are
%     a lot of NaNs. For example, interpolate 1:100 samples of LED2 positions.
%   * run fixLedAssignment.
%   * remove position jumps.
%
%  USAGE
%   fixedPos = general.fixLedAssignment(pos)
%   pos         Nx5 matrix of position samples, [t x y x2 y2].
%   fixedPos    Nx5 matrix of fixed position samples, [t x y x2 y2].
%
function fixedPos = fixLedAssignment(pos)
    costOfUnassignment = 100;
    pt1_lastGood = pos(1, 2:3);
    pt2_lastGood = pos(1, 4:5);

    fixedPos = pos;

    isDebug = false;
    if isDebug
        figure, axis([0 800 0 600]); grid on; hold on; %#ok<*UNRCH>
    end

    % process only points where at least one LED is not NaN. This should shorten
    % the loop, though will take time to find all these indices.
    validPoints = find(~isnan(pos(:, 2)) | ~isnan(pos(:, 4)));
    N = length(validPoints);
    
    w = 50; % Width of progress bar
    disp(['  0%[>', repmat(' ', 1, w), ']']);
    
    for i = 1:N
        idx = validPoints(i);
        pt1 = pos(idx, 2:3);
        pt2 = pos(idx, 4:5);

        points = [pt1; pt2];
        points(isnan(points)) = Inf;
        cost = zeros(2, 2);
        points_prev = [pt1_lastGood; pt2_lastGood];

        % create cost matrix directly without loops. Saves us around 1 sec.
        cost(1, 1) = pdist2(points(1, :), points_prev(1, :));
        cost(1, 2) = pdist2(points(2, :), points_prev(1, :));
        cost(2, 1) = pdist2(points(1, :), points_prev(2, :));
        cost(2, 2) = pdist2(points(2, :), points_prev(2, :));

        assignments = assignDetectionsToTracks(cost, costOfUnassignment);
        led1Idx = find(assignments(:, 1) == 1, 1);
        if isempty(led1Idx)
            pt1 = [nan nan];
        else
            pt1 = points(assignments(led1Idx, 2), :);
        end
        led2Idx = find(assignments(:, 1) == 2, 1);
        if isempty(led2Idx)
            pt2 = [nan nan];
        else
            pt2 = points(assignments(led2Idx, 2), :);
        end

        if ~isnan(pt1(1))
            pt1_lastGood = pt1;
        end
        if ~isnan(pt2(1))
            pt2_lastGood = pt2;
        end
        fixedPos(idx, 2:3) = pt1;
        fixedPos(idx, 4:5) = pt2;
        
        percent = 100 * i / N;
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);

        if isDebug
            plot(pt1(1), pt1(2), '.g');
            plot(pt2(1), pt2(2), '.r');
            title(sprintf('frame %u/%u', i, N));
            pause(0.03);
        end
    end
    disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
end
