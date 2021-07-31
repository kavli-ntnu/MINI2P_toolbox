% Interpolate position data
%
% Interpolates 1D or 2D position points on time interval. This function
% can be used to interpolate NaN values in position data.
%
%  USAGE
%   varargout = general.interpolatePositions(t, pos)
%   t            Vector of time measurements.
%   pos          One or two column matrix with position data.
%   varargout    Vectors of interpolated values each in separate variable.
%
function varargout = interpolatePositions(t, pos)
    if nargin < 2
        error('Incorrect number of parameters (type ''help <a href="matlab:help general.interpolatePositions">general.interpolatePositions</a>'' for details).');
    end

    if nargout < 1
        return;
    end
    
    warning('off', 'MATLAB:interp1:NaNstrip'); % disable warnings about NaN values in columns.
            % That's the point - we are interpolating NaN values.
    
    varargout{1} = interp1(t, pos(:, 1), t, 'pchip');
    if size(pos, 2) > 1 && nargout >= 2
        varargout{2} = interp1(t, pos(:, 2), t, 'pchip');
    end
    
    warning('on', 'MATLAB:interp1:NaNstrip');
end
