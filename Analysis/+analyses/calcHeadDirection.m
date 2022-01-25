% Calculate head direction.
%
% Calculates the head direction for each position sample pair. Direction
% is defined as east = 0 degrees, north = 90 degrees, west = 180 degrees,
% south = 270 degrees. Direction is set to NaN for missing samples.
% Position matrix contains information about front and back LED. Head
% direction is the counter-clockwise direction from back LED to the front.
%
%  USAGE
%    hd = calcHeadDirection(positions)
%
%    positions          Animal's position data, Nx5. Position data should
%                       contain timestamps (1 column), X/Y coordinates of
%                       first LED (2 and 3 columns correspondingly), X/Y
%                       coordinates of the second LED (4 and 5 columns
%                       correspondingly).
%                       it is assumed that positions(:, 2:3) correspond to
%                       front LED, and positions(:, 4:5) to the back LED.
%                       The resulting hd is the directon from back LED to
%                       the front LED.
%    hd                 Vector of head directions in degrees.
%
function hd = calcHeadDirection(positions)

    if size(positions, 2) < 5
        error('Position data should be 2D (type ''help <a href="matlab:help analyses.calcHeadDirection">analyses.calcHeadDirection</a>'' for details).');
    end
    % 1 column contains timestamps
    x1 = positions(:, 2);
    y1 = positions(:, 3);
    x2 = positions(:, 4);
    y2 = positions(:, 5);

    hd = rem(atan2d(y2-y1, x2-x1) + 180, 360);
end
