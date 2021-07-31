% Generate hexagonal grid
%
% This function generates a hexagonal grid pattern on a rectangle of specified dimension. 
%
%  USAGE
%   points = hexGrid(rect, size, plot)
%   rect    [x y w h] coordinates of a rectangle for the pattern. x, y - initial position,
%           w, h - width and height.
%   size    distance between hexes (default 1)
%   plot    True of false. If true, then the grid pattern is plotted on a new figure. Default is
%           false.
%   points  Points of the pattern. NB! Contains only unique points. There could be NaN values.
%

function points = hexGrid(rect, size, doPlot)
    
    [x, y, w, h] = deal(rect(1), rect(2), rect(3), rect(4));
    if nargin < 2
        size = 1;
        doPlot = false;
    end
    if nargin < 3
        doPlot = false;
    end
    [c, s] = deal(cos(pi/6), sin(pi/6));
    xlist = 0:size*(2+2*s):w-size*(1+s);
    ylist = 0:size*(c):h;
    offset_y = repmat(ylist, length(xlist), 1);
    offset_x = repmat(xlist, length(ylist), 1)';
    offset_x(:, 2:2:end) = offset_x(:, 2:2:end) + size*(1+s);
    offsets = exp(xlist') * exp(ylist);
    
    x0 = size * cumsum([0, s, 1, s, NaN]); % half hex
    y0 = size * [0 c c 0, NaN]; % half hex
    x3 = log(exp(x0') * exp(offset_x(:))');
    y3 = log(exp(y0') * exp(offset_y(:))');

    if doPlot
        figure;
        plot(x + x3(:), y + y3(:));
        axis equal; % optional
        hold on;
    end
    
    points(:, 1) = x + x3(:);
    points(:, 2) = y + y3(:);
    
    %% add centre fields
    y_step = ylist(2);
    numy = round(h / y_step);
    tmpy = (0:numy)*y_step;
    y_centres = tmpy(1:2:end);

    x_step = size/2;
    x_left = unique(offset_x) - x_step;
    x_centres = x_left(2:end);
    if length(xlist) > 1
        x_right = xlist(2) + x_left;
    else
        x_right = size*2 + x_step;
        x_right(2) = x_right(1);
    end
    x_centres = [x_centres; x_right(1:end-1)];
    x_centres = unique(x_centres);
    x_centres = x_centres(1:2:end);

    even_x_centres = x_centres + size + size/2;
    even_y_centres = y_centres + y_step;
    even_y_centres(end) = [];

    more_points = zeros(length(y_centres)*length(x_centres), 2);
    point_ind = 1;
    for i = 1:length(x_centres)
        for j = 1:length(y_centres)
            if doPlot
                plot(x_centres(i), y_centres(j), 'or');
            end
            more_points(point_ind, :) = [x_centres(i) y_centres(j)];
            point_ind = point_ind + 1;
        end
    end
    points = cat(1, points, more_points);
    
    more_points = zeros(length(even_x_centres) * length(even_y_centres), 2);
    point_ind = 1;
    for i = 1:length(even_x_centres)
        for j = 1:length(even_y_centres)
            if doPlot
                plot(even_x_centres(i), even_y_centres(j), 'or');
            end
            more_points(point_ind, :) = [even_x_centres(i) even_y_centres(j)];
            point_ind = point_ind + 1;
        end
    end
    points = cat(1, points, more_points);
    
    points = unique(points, 'rows');
end