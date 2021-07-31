% Calculate animals' speed based on position samples
%
% This script calculates momentary speed and smooths it. The result shall be very similar
% to the result of using a Kalman filter. However, this function is much faster.
% Lowess method is used for smoothing with span of 0.8 seconds.
%
%  USAGE
%   v = general.speed(pos)
%   pos         Matrix with position data. Should be at least of size Nx2 and
%               structure [t x]. More general structure is [t x y].
%   v           Vector/matrix of size N or NxM with speed values. Speed units depends on units
%               of position samples. If positions are in cm, then [v] is [cm/s].
%               If pos contains multiple positions groups (LEDs), then M is the number of groups.
%               Each column of v correponds to a position group, i.e. v(:, 1) is first LED,
%               v(:, 2) is speed for the second LED. Maximum number of two LEDs is supported.
%  SEE
%
%    See also smooth.
%
function v = speed(pos)
    if nargin < 1
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help general.speed">general.speed</a>'' for details).');
    end

    if size(pos, 2) < 2
        error('BNT:arg', 'Incorrect argument ''pos'' (type ''help <a href="matlab:help general.speed">general.speed</a>'' for details).');
    end

    diffs = diff(pos);
    dt = diffs(:, bntConstants.PosT);
    dx = diffs(:, bntConstants.PosX);
    if size(pos, 2) > 2
        dy = diffs(:, bntConstants.PosY);
    else
        dy = zeros(size(pos, 1)-1, 1);
    end

    Fq = 1/mean(dt);
    filterSpan = 0.8 * Fq;
    if size(pos, 2) > 3
        numLeds = 2;
    else
        numLeds = 1;
    end
    N = size(pos, 1);
    v = zeros(N, numLeds);
    
    if N <= filterSpan
        % output unsmoothed speed
        t = pos(:, bntConstants.PosT);
        x = pos(:, bntConstants.PosX);
        if size(pos, 2) > 2
            y = pos(:, bntConstants.PosY);
        else
            y = ones(1, size(pos, 1));
        end
        for i = 2:N-1
            v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
        end
        v(1, 1) = v(2, 1);
        v(end, 1) = v(end-1, 1);
        
        if numLeds > 1
            x = pos(:, bntConstants.PosX2);
            y = pos(:, bntConstants.PosY2);
            for i = 2:N-1
                v(i, 2) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
            end
            v(1, 2) = v(2, 2);
            v(end, 2) = v(end-1, 2);
        end
    else
        v_x = smooth(dx ./ dt, filterSpan, 'lowess');
        v_y = smooth(dy ./ dt, filterSpan, 'lowess');
        tmpV = sqrt(v_x.^2 + v_y.^2);
        tmpV(end+1) = tmpV(end);
        v(:, 1) = tmpV;
        if numLeds > 1
            dx = diffs(:, bntConstants.PosX2);
            dy = diffs(:, bntConstants.PosY2);
            v_x = smooth(dx ./ dt, filterSpan, 'lowess');
            v_y = smooth(dy ./ dt, filterSpan, 'lowess');
            tmpV = sqrt(v_x.^2 + v_y.^2);
            tmpV(end+1) = tmpV(end);
            v(:, 2) = tmpV;
        end
    end
end
