% Extrapolates place fields on a map
%
% This function extrapolates place fields of a map in a vertical direction. The extrapolation
% is done based on grid properties and is carried out by means of creating parallel axes
% of a grid cell. The intersection of axes define new place fields.
%
%  USAGE
%   extendedMap = analyses.extendMapByAxes(map, fields, gridStat)
%   map             A rate map structure (output of analyses.map)
%   fields          Detected place fields (output of analyses.placefield)
%   gridStat        Grid statistics (output of analyses.gridnessScore)
%   extendedMap     A 2D matrix
%

function extendedMap = extendMapByAxes(map, fields, gridStat)
    if isempty(fields)
        extendedMap = [];
        return;
    end
    vert = [[fields(:).peakX]' [fields(:).peakY]']; % assume that vert is already sorted by x

    mapSize = size(map.z);
    minSize = min(size(map.z));
    merged = zeros(minSize*2);
    if mapSize(1) < mapSize(2)
        merged(minSize+1:minSize*2, :) = map.z;
    else
        error('BNT:err:dimensions', 'Only rectangular maps with width longer than height are supported!');
    end

    fieldRadii = sqrt([fields(:).area] / pi);
    meanFieldRadius = median(fieldRadii) * 2; % 2 is an arbitrary number. This value is used to create new place fields.
                                              % Higher radius value results in better fields with helpers.gauss2d
    averagePeak = median([fields(:).peak]) * 1.2; % 1.2 is an arbitrary number. Higher peak frequency results in better fields.

    % DEBUG: view intermediate results
    % h = figure(); plot.colorMap(merged), axis image, hold on;
    % xlim([-20 mapSize(2) + 20]); % increase the limits to see axis lines beyond the original map
    % ylim([-20 mapSize(2) + 20]);

    allLines = [0 0 0 0]; % all created axis lines go here
    k = 1; % index for allLines

    % ====================
    % We have three different orientation values from gridStat. Axes (parallel lines) are created based on these values.
    % However, only two orientation values are used. If third is used, then I found that lines do not intersect in a single
    % point. Thus additional processing is needed to detect close points and to somehow merge them.
    % ====================

    vertInd = 1; % we create the axes relative to (starting from) this place field
    newPointRho = 10; % 10 is an arbitrary radius. Could be any number

    %% process the first orientation
    initPoint = vert(vertInd, :) + [0 minSize]; % start from the very first field, add minSize because the point is
                                          % on the merged map now
    theta = gridStat.orientation(1) + 90; % orthogonal direction
    stepX = gridStat.spacing(2)*cosd(theta); % one step in x direction
    stepY = gridStat.spacing(2)*sind(theta); % one step in y direction
    allX = initPoint(1) - stepX:stepX:size(merged, 2) + stepX; % create even more points than we need
                                                               % to cover areas outside the map
    allY = initPoint(2) + stepY:-stepY:0 - stepY;
    if length(allY) > length(allX)
        allY = allY(1:length(allX));
    end
    if length(allX) > length(allY)
        allX = allX(1:length(allY));
    end

    for i = 1:length(allX)
        initPoint = [allX(i) allY(i)];
        theta = gridStat.orientation(1);
        x = initPoint(1) + newPointRho*cosd(theta);
        y = initPoint(2) + newPointRho*sind(theta);
        line = createLine(initPoint, [x y]);
        allLines(k, :) = line;
        k = k + 1;

        % DEBUG: view intermediate results
        % drawLine(line, 'color', 'w', 'linewidth', 2);
    end

    if gridStat.orientation(3) < 90
        %% process the third orientation
        initPoint = vert(vertInd, :) + [0 minSize];
        theta = gridStat.orientation(3) - 90; % orthogonal direction
        stepX = gridStat.spacing(2)*cosd(theta); % one step in x direction
        stepY = gridStat.spacing(2)*sind(theta); % one step in y direction
        allX = initPoint(1) - stepX:stepX:size(merged, 2) + stepX; % create even more points than we need
                                                                   % to cover areas outside the map
        allY = initPoint(2) + abs(stepY):stepY:0 + stepY;
        if length(allX) > length(allY)
            allX = allX(1:length(allY));
        end
        if length(allY) > length(allX)
            allY = allY(1:length(allX));
        end

        for i = 1:length(allX)
            initPoint = [allX(i) allY(i)];
            theta = gridStat.orientation(3);
            x = initPoint(1) + newPointRho*cosd(theta);
            y = initPoint(2) + newPointRho*sind(theta);
            line = createLine(initPoint, [x y]);
            allLines(k, :) = line;
            k = k + 1;

            % DEBUG: view intermediate results
            % drawLine(line, 'color', 'r', 'linewidth', 2);
        end
    else
        % try to process the second orientation, because when gridStat.orientation(3) > 90,
        % then the resulting lines are parallel to the first orientation
        initPoint = vert(vertInd, :) + [0 minSize];
        theta = gridStat.orientation(2) + 90; % orthogonal direction
        stepX = gridStat.spacing(2) * cosd(theta);
        stepY = gridStat.spacing(2) * sind(theta);

        factorX = abs(ceil(initPoint(1)/stepX));
        factorY = abs(ceil(initPoint(2)/stepY));

        allX = initPoint(1) + stepX*factorX:-stepX:mapSize(2)-stepX;
        allY = initPoint(2) - stepY*factorY:stepY:mapSize(2)+stepY;
        if length(allX) > length(allY)
            allX = allX(1:length(allY));
        end
        if length(allY) > length(allX)
            allY = allY(1:length(allX));
        end

        for i = 1:length(allX)
            initPoint = [allX(i) allY(i)];
            theta = gridStat.orientation(2);
            x = initPoint(1) + newPointRho*cosd(theta);
            y = initPoint(2) + newPointRho*sind(theta);
            line = createLine(initPoint, [x y]);
            allLines(k, :) = line;
            k = k + 1;

            % DEBUG: view intermediate results
            % drawLine(line, 'color', 'y', 'linewidth', 2);
        end
    end

    %% detect intersections
    mu = zeros(1, 2);
    amplitudes = 0;
    sigmas = zeros(1, 2);

    k = 1;
    for i = 1:size(allLines, 1)
        line1 = allLines(i, :);
        for j = i+1:size(allLines, 1)
            line2 = allLines(j, :);
            pt = intersectLines(line1, line2, 1e-10);
            if ~any(isinf(pt)) && ~any(isnan(pt))
                % a valid point
                mu(k, :) = pt;
                if pt(1) > size(merged, 2) || pt(2) > size(merged, 1) || pt(1) < 0 || pt(2) < 0
                    % make fields that are outside the map even bigger, so that their border is more likely visible
                    sigmas(k) = meanFieldRadius * 1.5; % 1.5 is an arbitary number
                    amplitudes(k, :) = averagePeak * 1.5; %#ok<AGROW>
                else
                    sigmas(k) = meanFieldRadius;
                    amplitudes(k, :) = averagePeak; %#ok<AGROW>
                end

                k = k + 1;
            end
        end
    end
    % DEBUG: view intermediate results
    % close(h);

    %% create the resulting map
    extendedMap = helpers.gauss2d(mu, sigmas', size(merged), amplitudes);
end

