% Calculate animal angular head velocities
%
% Angular head velocity is calculated according to H.T. Blair and P. E. Sharp:1995 paper.
% The momentary angular velocity of the animal's head is calculated as the difference
% in the angle of head direction between successive samples.
%
%  USAGE
%   v = analyses.angularHeadVelocity(headDirections)
%   headDirections      matrix Nx2 of animal's head directions in degrees.
%                       1 column contains  timestamps.
%                       2 column contains head directions. See analyses.calcHeadDirection.
%                       N must be greater or equal to 2.
%   v                   vector of calculated angular head velocities. [degrees/sec] * 1/<sample_time>
%
function v = angularHeadVelocity(headDirections)
    t = headDirections(:, 1);
    if length(t) < 2
        error('BNT:args:length', 'headDirections must have at least 2 samples, i.e. N >= 2.');
    end
    dt = diff(t);
    t_diff = t(1:end-1) + dt/2;

    % take derivative
    d0 = diff(headDirections(:, 2)) ./ dt;

    % interpolate it to the whole time range
    d1 = interp1(t_diff, d0, t(2:end, 1));

    % add last interpolated sample, so that size(v) = size(headDirections(:, 2))
    v = [d0; d1(end)]; % v units are [degrees/sec]. And actual value
            % can be seen as v * <sample_time_coefficient> [deg/sec].

    v = v * -1; % negative turns are clockwise
end
