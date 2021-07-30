% Map z on (x, y) where x, y and z are time-varying variables (samples)
%
% This function computes a map where one time-varying variable z is represented as a function
% of one or two time-varying variables x and y. The variable z is typically a list of spike
% timestamps, and x and y are spatial coordinates.
%
%  USAGE
%    map = analyses.map(pos, z, <options>)
%    pos            Matrix with positions samples in form [t x y] or [t x], where
%                   t              timestamps for x and y
%                   x              x values
%                   y              optional y values
%                   If pos contains more than 3 columns, then only first 3 are used.
%    z              list of timestamps or Nx2 matrix of time-depended variables. One example
%                   could be Nx2 matrix of instanteneous animal speed in form [t s].
%                   2 dimensional z is only supported for 2D pos!
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 1).
%     'binWidth'    width of horizontal and vertical bins (default = [2.5]). If only one value is provided
%                   then the same bin width is used in both directions. The units of binWidth is the same
%                   as the units of positions. If you have your positions in cm, then binWidth is in cm.
%     'limits'      Vector with 2 or 4 elements. Minimum and maximum values of x and y. Default value
%                   is [nanmin(x) nanmax(x) nanmin(y) nanmax(y)]. This value is used to arrange
%                   data in bins. For example, if you want to calculate two firing rate maps with
%                   common bin-structure you should use identical limits parameter.
%     'minTime'     minimum time spent in each bin (in s, default = 0).
%     'maxGap'      z values recorded during time gaps between successive (x,y) samples exceeding
%                   this threshold (e.g. undetects) will not be interpolated; also, such long gaps
%                   in (x,y) sampling will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s).
%     'type'        'linear' if z is linear (default), 'circular' otherwise.
%     'posGroups'   Position group indices. Positions can be separated in groups. This value should
%                   be a vector of size length(t). posGroups(i) == n, means that position sample i
%                   belongs to group number n. Default value all ones, meaning that positions belong
%                   to single group. See groupType.
%     'spkGroups'   Spike group indices. Spikes can be separated in groups. This value should be
%                   a vector of size length(z). spkGroups(i) == n, means that spike timestamp i
%                   belongs to group number n. Default value all ones, meaning that positions belong
%                   to single group. See groupType.
%     'groupType'   'joint' if posGroups and spkGroups were grouped in similar manner (for example,
%                   grouped by 1 sec intervals), 'separate' otherwise. You can get separate groups
%                   if you process several cells simultaneously. That's it you have  one position
%                   group and several spike groups. Default is 'joint'.
%     'blanks'      'on' results on NaN values in map for unvisited places (regardless of smoothing).
%                   'off' unvisited places will have a very tiny firing rate. Default is 'on'
%    =========================================================================
%
%  OUTPUT
%
%    Array of structures:
%    map.x          x bins
%    map.y          y bins
%    map.z          rate map
%    map.count      count map, for 1 dimensional z this is the count of events (spikes) happened in each bin.
%                   For 2 dimensional z, this is the number of times animal has been in each bin.
%    map.time       occupancy map (in s), how much time animal has spent in each bin.
%    map.zRaw       If smooth is > 0, returns z without smoothing
%    map.countRaw   Raw (unsmoothed) version of map.count
%    map.timeRaw    Raw (unsmoothed) version of map.time
%    map.zBins      Linear index of z bin belonging. In general, contains bin number to which every z-sample
%                   belongs.
%
%    Number of elements in output array depends on position and spike groups. If groups are separate,
%    then resulting array will have <number of spike groups> * <number of position groups> elements.
%    If they are joint, then output array will have <number of position groups> elements.
%
%  NOTES
%
%    x values are arranged in columns and y values in rows in all output matrices (e.g. 'map.z').
%
%  SEE
%
%    See also plot.colorMap, general.accumulate.

% Copyright (C) 2002-2011 by Michaël Zugaro, 2013-2015 Vadim Frolov
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function map = map(v, z, varargin)
    % Check parameter sizes
    if size(v, 2) < 2
        error('Parameter ''pos'' should have at least 2 columns (type ''help <a href="matlab:help analyses.map">analyses.map</a>'' for details).');
    end
    if (size(z, 2) < 1 || size(z, 2) > 2) && ~isempty(z)
        error('Parameter ''z'' should have 1 or 2 columns (type ''help <a href="matlab:help analyses.map">analyses.map</a>'' for details).');
    end
    if size(z, 2) == 2 && ~isequal(size(z, 1), size(v, 1))
        error('Parameter ''z'' should have the same number of rows as ''pos''. (type ''help <a href="matlab:help analyses.map">analyses.map</a>'' for details).');
    end

    inp = inputParser;

    % Default values
    defaultMaxGap = 0.1;
    defaultSmooth = 1;
    defaultBinWidth = [2.5 2.5];
    defaultMinTime = 0;
    defaultType = 'linear';
    defaultPosIndices = ones(size(v, 1), 1);
    defaultSpkIndices = ones(size(z, 1), 1);
    defaultGroupType = true;
    defaultLimits = [];

    checkSmooth = @(x) helpers.isdvector(x, '>=0') && length(x) <= 2;
    checkScalar = @(x) helpers.isdscalar(x, '>=0');
    checkType = @(x) helpers.isstring(x, 'circular', 'linear');
    checkVector = @(x) helpers.isivector(x, '>0');
    checkGroupType = @(x) helpers.isstring(x, 'joint', 'separate');
    checkBlanks = @(x) helpers.isstring(x, 'on', 'off');
    checkJustVector = @(x) isvector(x);

    addRequired(inp, 'v');
    addRequired(inp, 'z');
    addParameter(inp, 'smooth', defaultSmooth, checkSmooth);
    addParameter(inp, 'binWidth', defaultBinWidth, checkSmooth);
    addParameter(inp, 'minTime', defaultMinTime, checkScalar);
    addParameter(inp, 'maxGap', defaultMaxGap, checkScalar);
    addParameter(inp, 'type', defaultType, checkType);
    addParameter(inp, 'posGroups', defaultPosIndices, checkVector);
    addParameter(inp, 'spkGroups', defaultSpkIndices, checkVector);
    addParameter(inp, 'groupType', defaultGroupType, checkGroupType);
    addParameter(inp, 'blanks', 'on', checkBlanks);
    addParameter(inp, 'limits', defaultLimits, checkJustVector);

    parse(inp, v, z, varargin{:});
    % get parsed results
    smooth = inp.Results.smooth;
    binWidth = inp.Results.binWidth;
    minTime = inp.Results.minTime;
    maxGap = inp.Results.maxGap;
    type = inp.Results.type;
    posIndices = inp.Results.posGroups;
    spkIndices = inp.Results.spkGroups;
    areGroupsJoint = inp.Results.groupType;
    sBlanks = inp.Results.blanks;
    showBlanks = strcmpi(sBlanks, 'on');
    limits = inp.Results.limits;

    %% rest of the function
    if isempty(v)
        return;
    end

    % Some info about x, y and z
    pointProcess = (isempty(z) | size(z, 2) == 1);
    t = v(:,1);
    x = v(:,2);
    if size(v, 2) >= 3
        y = v(:, 3);
    else
        y = [];
    end

    if ~isempty(limits) && ~isempty(y) && length(limits) < 4
        % in 2D case and if limits are provided, there should be 4 elements
        error('Incorrect value for property ''limits'' (type ''help <a href="matlab:help analyses.map">analyses.map</a>'' for details).');
    end

    % Number of bins for x and y
    binWidthX = binWidth(1);
    if length(binWidth) == 1
        binWidthY = binWidthX;
        binWidth(2) = binWidth; %#ok<NASGU>
    else
        binWidthY = binWidth(2);
    end

    originalX = x;
    originalY = y;
    originalT = t;
    originalZ = z;

    numPosGroups = length(unique(posIndices));
    numSpikeGroups = length(unique(spkIndices));
    if areGroupsJoint
        numMaps = numPosGroups;
    else
        numMaps = numPosGroups * numSpikeGroups;
    end

    map(numMaps).y = [];
    map(numMaps).x = [];
    map(numMaps).count = [];
    map(numMaps).time = [];
    map(numMaps).z = [];
    map(numMaps).Nspikes = [];
    map(numMaps).zRaw = [];
    map(numMaps).countRaw = [];
    map(numMaps).timeRaw = [];
    map(numMaps).peakRate = [];
    map(numMaps).minRate = [];
    if isempty(z)
        return
    end

    m = 1; % linear index of current map
    for p = 1:numPosGroups
        x = originalX(posIndices == p);
        if ~isempty(originalY)
            y = originalY(posIndices == p);
        end
        t = originalT(posIndices == p);
        if areGroupsJoint
            z = originalZ(spkIndices == m, :);
        end

        if isempty(limits)
            limitsX = [nanmin(x) nanmax(x)]; % TODO: get limits from the calibration
        else
            limitsX = [limits(1) limits(2)];
        end
        nBinsX = ceil((limitsX(2) - limitsX(1)) / binWidthX);
        if nBinsX == 0
            continue;
        end

        if ~isempty(y)
            if isempty(limits)
                limitsY = [nanmin(y) nanmax(y)]; % TODO: get limits from the calibration
            else
                limitsY = [limits(3) limits(4)];
            end
            nBinsY = ceil((limitsY(2) - limitsY(1)) / binWidthY);
            if nBinsY == 0
                continue;
            end
        end

        % Bin x and y
        edges = limitsX(1):binWidthX:limitsX(1) + binWidthX*nBinsX; % create nBinsX bins
        [~, ~, x] = histcounts(x, edges);
        x(x == 0) = nan; % x == 0 indicates values outside edges range, could happen depending on the limitsX

        if ~isempty(y)
            edges = limitsY(1):binWidthY:limitsY(1) + binWidthY*nBinsY;
            [~, ~, y] = histcounts(y, edges);
            y(y == 0) = nan;
        end

        % Duration for each (X,Y) sample (clipped to maxGap)
        dt = diff(t);
        if isempty(dt)
            dt = 0;
        end
        dt(end+1) = dt(end); %#ok<AGROW>
        dt(dt>maxGap) = maxGap;

        xSpace = limitsX(1):binWidthX:limitsX(1) + binWidthX*nBinsX;
        if isempty(y)
            time = general.accumulate(x, dt, nBinsX)'; % how much time animal spent at a specific bin
            time = round(time * 10e6) / 10e6; % remove floating numbers artefacts
        else
            ySpace = limitsY(1):binWidthY:limitsY(1) + binWidthY*nBinsY;
            % accumulate returns matrix MxN where M are columns and not rows, as opposite to
            % standard Matlab matrix notation.
            time = general.accumulate([x y], dt, [nBinsX nBinsY])'; % how much time animal spent at a specific bin
            time = round(time * 10e6) / 10e6; % remove floating numbers artefacts
        end

        if areGroupsJoint
            numIter = 1;
        else
            numIter = numSpikeGroups;
        end

        for s = 1:numIter
            if isempty(time) || isempty(z)
                % sometimes time can be empty if all position sampels are NaN
                if isempty(y)
                    map(m).z = zeros(1, nBinsX);
                else
                    map(m).z = zeros(nBinsY, nBinsX);
                end

                m = m + 1; % increase map linear index

                continue;
            end

            if ~areGroupsJoint
                z = originalZ(spkIndices == s, :);
            end

            if pointProcess
                % Count occurrences for each (x,y) timestamp, i.e. spike spread across positions
                n = histc(z, t);
            else
                % Interpolate z at (x,y) timestamps
                %[z, discarded] = Interpolate(z, t, 'maxGap', maxGap);
                %if isempty(z), return; end
                if strcmp(type, 'circular')
                    range = helpers.isradians(z(:, 2));
                    z(:, 2) = exp(1i*z(:, 2));
                end
                n = 1;
            end

            % Computations
            if isempty(y)
                % 1D (only x)
                map(m).x = xSpace;

                if isempty(z)
                    Nspikes = zeros(nBinsX, 1);
                else
                    Nspikes = general.accumulate(x, n, nBinsX);
                end

                map(m).countRaw = Nspikes'; % how many times (# spikes) animal was at a specific bin
                map(m).timeRaw = time;
                map(m).Nspikes = Nspikes'; % used elsewhere

                if pointProcess
                    zRaw = map(m).countRaw ./ (map(m).timeRaw + eps);
                    map(m).zRaw = zRaw;
                else
                    % zRaw contains sum of z(:, 2) elements that fall in a bin
                    zRaw = general.accumulate(x, z(:, 2), nBinsX)';
                    % divide sum of elements by number of visits to get mean value
                    zRaw = zRaw ./ (map.countRaw + eps);
                    map(m).zRaw = zRaw;
                end

                map(m).z = general.smooth(zRaw, smooth);
                map(m).time = general.smooth(map(m).timeRaw, smooth);
                map(m).count = general.smooth(map(m).countRaw, smooth);

                if ~isequal(size(map(m).z), size(zRaw))
                    map(m).z = map(m).z';
                    map(m).time = map(m).time';
                    map(m).count = map(m).count';
                end
            else
                % 2D (x and y)
                map(m).x = xSpace;
                map(m).y = ySpace;

                % Nspikes = general.accumulate([x y], n, nBins);
                % map(m).Nocc = general.smooth(general.accumulate([x y]), smooth);

                if isempty(z)
                    Nspikes = zeros([nBinsX nBinsY]);
                else
                    % for non-point process this is number of animal visits in a bin
                    Nspikes = general.accumulate([x y], n, [nBinsX nBinsY]);
                end
                map(m).countRaw = Nspikes'; % how many times (# spikes) animal was at a specific bin
                map(m).timeRaw = time;
                map(m).Nspikes = Nspikes';

                if pointProcess
                    zRaw = map(m).countRaw ./ (map(m).timeRaw + eps);
                    map(m).zRaw = zRaw;
                else
                    % zRaw contains sum of z(:, 2) elements that fall in a bin
                    zRaw = general.accumulate([x y], z(:, 2), [nBinsX nBinsY])';
                    % divide sum of elements by number of visits to get mean value
                    zRaw = zRaw ./ (map.countRaw + eps);
                    map(m).zRaw = zRaw;
                end

                map(m).z = general.smooth(zRaw, smooth);
                map(m).time = general.smooth(map(m).timeRaw, smooth);
                map(m).count = general.smooth(map(m).countRaw, smooth);
            end

            if showBlanks
                selected = time == 0;
                map(m).time(selected) = NaN;
                map(m).z(selected) = NaN;
                map(m).zRaw(selected) = NaN;
                map(m).countRaw(selected) = NaN;
                map(m).timeRaw(selected) = NaN;
            end

            % Circular z
            if strcmp(type, 'circular')
                map(m).z = general.wrap(angle(map(m).z), range);
            end

            % Discard regions with insufficient sampling
            ind = map(m).time < minTime;
            if showBlanks
                map(m).z(ind) = NaN;
                map(m).zRaw(ind) = NaN;
            else
                map(m).z(ind) = 0;
                map(m).zRaw(ind) = 0;
            end
            map(m).peakRate = nanmax(nanmax(map(m).z));
            map(m).minRate = nanmin(nanmin(map(m).z));
            m = m + 1;
        end
    end
end % parent function

