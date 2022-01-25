% Locate place fields in a 1D firing map.
%
% Locates the place fields and calculates the field sizes and spacing
% between fields.
%
%  USAGE
%   [fieldsCurve, fields, spacing] = analyses.placefield1D(map, <options>)
%   map           Firing rate map either structure obtained using <a href="matlab:help analyses.map">analyses.map</a>
%                 or a matrix that represents firing map. Note that if you want to filter place fields based on number
%                 of spikes, then you should use map structure.
%   <options>     optional list of property-value pairs (see table below)
%
%   ======================================================================================
%    Properties     Values
%   --------------------------------------------------------------------------------------
%    'thresholdType' Type of the firing rate threshold. '%' specifies that 'threshold'
%                    value is in range [0; 1] and rates lower than 'threshold * curveAmplitude'
%                    are not considered as place fields. curveAmplitude is the range of a firing
%                    curve, i.e. max(curve) - min(curve).
%                    'rate' specifies that 'threshold' value is an absolute rate. Values lower than
%                    'threshold' are not considered as place fields. Default it '%'.
%
%    'threshold'    Threshold for a firing rate. How the value is treated depends on options
%                   'thresholdType'. Default = 0.2 and 'thresholdType' = '%'.
%
%    'minRows'      Minimum number of row bins in a place field. Fields with
%                   fewer bins are not considered as place fields. Remember to
%                   adjust this value when you change the bin width. Default = 3.
%                   If minRows == 0, then minimum number of bins is not checked.
%
%    'binWidth'     Bin length in centimetres. Used to calculate size of fields in cm.
%                   Default is 1.
%
%    'minSpikes'    Minimum number of spikes in a place field. Default is 0, which means
%                   that minimum number of spikes is not checked. To use this value 'map'
%                   must be a result of function <a href="matlab:help analyses.map">analyses.map</a>.
%
%    'minDistance'  Minimum number of bins between adjacent fields. Distance is calculated between
%                   end of one field and start of another field. If distance is smaller or
%                   equals to 'minDistance', then two fields are merged together. Note that more
%                   than two fields can be merged if they are all adjacent. Default is 0, which
%                   means that fields should have a common border.
%
%    'debug'        'on' plot turning curve, place fields and their centres of mass. 'off' do not
%                   plot anything. Default is 'off'.
%
%    'pos'          Position samples. Used to calculate posInd. If not provided,
%                   then posInd will be an empty matrix.
%   ======================================================================================
%
%   fieldsCurve     Vector of the same size as underline firing rate curve. The elements of fieldsCurve
%                   are integer values greater or equal to 0. 0 elements correspond to background (not
%                   a place field). Values of 1 constitute first field. Values of 2 constitute second
%                   field; and so on.
%   fields          Structure with information about each field. Structure fields are:
%       col         vector of columns that constitute this field.
%       row         vector of rows (artificial). Same size as col, all ones;
%       peak        field peak;
%       x, y        field centre of mass point. See NOTES for details about centre of mass calculation;
%       meanRate    mean firing rate;
%       width       field width in cm;
%       size        size of the field. Do not rely on it! It's from old code and could be misleading.
%                   Calculated as width * <mean firing rate of a field>.
%       posInd      Indices of position samples that correspond to this field. In case position sample
%                   matrix is of size Nx5 ([t x y x1 y1]), posInd corresponds to the left most position
%                   columns ([x y]). If pos argument is not provided, then posInd will be an empty matrix.
%   spacing         Matrix of spacing between fields. Entry (i, j) equals to distance between fields i and j.
%
%  NOTE
%   Since centre of mass for a curve is somewhat undefined, it is calculated for plate bounded by a firing curve
%   and a function y = 0 for all x. See http://tutorial.math.lamar.edu/Classes/CalcII/CenterOfMass.aspx
%
function [fieldsCurve, fields, spacing] = placefield1D(mapS, varargin)
    if nargin < 1 || mod(length(varargin), 2) ~= 0
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
    end

    % Default values
    threshold = 0.1;
    minRows = 2;
    binWidth = 1; % [cm];
    minSpikes = 0;
    minFieldDistance = 0; % [bins]
    debug = false;
    isThrshPercentage = true;
    pos = [];

    % Parse parameter list
    i = 1;
    while i < length(varargin)
        if ~ischar(varargin{i}),
            error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).']);
        end

        switch(lower(varargin{i})),
            case 'threshold'
                threshold = varargin{i+1};
                i = i + 2;

            case 'minrows'
                minRows = varargin{i+1};
                if ~helpers.isdscalar(minRows, '>=0')
                    error('Incorrect value for property ''minRows'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            case 'binwidth'
                binWidth = varargin{i+1};
                if ~helpers.isdscalar(binWidth, '>0')
                    error('Incorrect value for property ''binWidth'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            case 'minspikes'
                minSpikes = varargin{i+1};
                if ~helpers.isdscalar(minSpikes, '>=0')
                    error('Incorrect value for property ''minSpikes'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                if minSpikes > 0 && ~isstruct(mapS)
                    error('You should only use minSpikes with a structure-like firing curve (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            case 'mindistance'
                minFieldDistance = varargin{i+1};
                if ~helpers.isdscalar(minFieldDistance, '>=0')
                    error('Incorrect value for property ''minDistance'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            case 'debug'
                debug = strcmpi(varargin{i+1}, 'on');
                i = i + 2;

            case 'thresholdtype'
                isThrshPercentage = strcmpi(varargin{i+1}, '%');
                if ~helpers.isstring(varargin{i+1}, '%', 'rate'),
                    error('Incorrect value for property ''thresholdType'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            case 'pos'
                pos = varargin{i+1};
                if ~ismatrix(pos) || size(pos, 2) < 2
                    error('Incorrect value for property ''pos'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
                end
                i = i + 2;

            otherwise,
                error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).']);
        end
    end
    if isThrshPercentage
        if ~helpers.isdscalar(threshold, '>=0', '<=1')
            error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
        end
    else
        if ~helpers.isdscalar(threshold, '>=0')
            error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
        end
    end

    if isstruct(mapS)
%         mapAxis = mapS.x;
        if size(mapS.z, 1) > 1 && size(mapS.z, 2) > 1
            error('Incorrect value for property ''map''. Firing curve is 2D, but should be a vector (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
        end

        curve = mapS.z;
    else
        if size(mapS, 1) > 1 && size(mapS, 2) > 1
            error('Incorrect value for property ''map''. It should be a vector (type ''help <a href="matlab:help analyses.placefield1D">analyses.placefield1D</a>'' for details).');
        end
%         mapAxis = 1:length(mapS);

        curve = mapS;
    end
    curveLen = length(curve);

    fieldsCurve = zeros(1, curveLen);
    fields = struct('row', {}, 'col', {}, ...
        'size', {}, 'peak', {}, ...
        'x', {}, 'y', {}, ...
        'meanRate', {}, 'width', {}, ...
        'posInd', {} ...
        );
    spacing = [];

    if curveLen == 0
        return
    end
    
    % extend curve, so that findpeaks function finds all peaks
    map = zeros(1, curveLen + 2);
    map(2:2+curveLen-1) = curve;

    [~, locs] = findpeaks(map);
    dMap = diff(map);
    dMap(end+1) = dMap(end);

    if isThrshPercentage
        threshold = min(curve) + ((max(curve) - min(curve)) * threshold);
    end

    startField = zeros(1, 100);
    stopField = zeros(1, 100);
    fieldWidth = zeros(1, 100);
    fieldSize = zeros(1, 100);
    curField = 1;

    for i = 1:length(locs)
        if map(locs(i)) < threshold
            continue;
        end
        curStopField = locs(i);
        while curStopField < length(dMap) && dMap(curStopField) < 0 && map(curStopField) >= threshold
            curStopField = curStopField + 1;
        end
        % substract one because of the loop
        curStopField = curStopField -1;

        curStartField = locs(i) - 1;
        while curStartField > 0 && dMap(curStartField) > 0 && map(curStartField) >= threshold
            curStartField = curStartField - 1;
        end
        if curStartField == 0
            curStartField = 1;
        end
        curStopField = curStopField - 1; % substract one cause of map -> curve translation

        if curStopField > length(curve)
            curStopField = length(curve);
        end

        startField(curField) = curStartField;
        stopField(curField) = curStopField;

        fieldWidth(curField) = (curStopField - curStartField + 1) * binWidth;
        fieldSize(curField) = fieldWidth(curField) * nanmean(curve(curStartField:curStopField));
        curField = curField + 1;
    end

    startField(curField:end) = [];
    stopField(curField:end) = [];
    fieldWidth(curField:end) = [];
    fieldSize(curField:end) = [];

    nFields = length(fieldWidth);
    removeField = false(1, nFields);

    fieldsCurve = zeros(1, curveLen);
    if nFields == 0
        fields = struct();
        fields(1) = [];
        spacing = [];
        return;
    end

    % merge neighbouring fields
    curField = 1;
    for i = 1:nFields
        if i + 1 > nFields
            break;
        end

        fieldsDist = abs(startField(i+1) - stopField(i));
        if fieldsDist <= minFieldDistance
            stopField(curField) = stopField(i+1);
            removeField(i+1) = true;
        else
            curField = i + 1;
        end
    end
    startField(removeField) = [];
    stopField(removeField) = [];
    fieldWidth(removeField) = [];
    fieldSize(removeField) = [];

    if ~isempty(find(removeField, 1))
        % adjust fields width and size after merge
        nFields = length(startField);
        for i = 1:nFields
            fieldWidth(i) = (stopField(i) - startField(i) + 1) * binWidth;
            fieldSize(i) = fieldWidth(i) * nanmean(curve(startField(i):stopField(i)));
        end
    end

%     fields(nFields) = struct('row', {}, 'col', {}, ...
%         'size', {}, 'peak', {}, ...
%         'x', {}, 'y', {}, ...
%         'meanRate', {}, 'width', {}, ...
%         'posInd', {} ...
%         );
%     fields(nFields) = struct('row', [], 'peak', 0, 'x', -1, 'y', -1, 'meanRate', -1, 'posInd', []);

    curField = 1;
    removeField = false(1, nFields);
    for i = 1:nFields
        fieldData = curve(startField(i):stopField(i));
        peak = max(fieldData);

        % check field peak rate
        if peak < threshold
            removeField(i) = true;
            continue;
        end

        % Number of spikes in the field
        if minSpikes > 0
            numSpikes = sum(mapS.Nspikes(startField(i):stopField(i)));
            if numSpikes < minSpikes
                % Marks field as having to few spikes
                removeField(i) = true;
                continue;
            end
        end

        % check number of bins
        if length(startField(i):stopField(i)) <= minRows
            removeField(i) = true;
            continue;
        end

        fields(curField).col = startField(i):stopField(i);
        fields(curField).row = ones(1, length(fields(curField).col));
        fields(curField).peak = peak;
        fields(curField).meanRate = nanmean(fieldData);
        fields(curField).width = fieldWidth(i);
        fields(curField).size = fieldSize(i);
        fieldsCurve(startField(i):stopField(i)) = curField;

        if ~isempty(pos)
            if isstruct(mapS)
                xSpace = mapS.x;
            else
                nBins = length(curve);
                limitsX = [nanmin(pos(:, bntConstants.PosX)) nanmax(pos(:, bntConstants.PosX))];
                xSpace = linspace(limitsX(1), limitsX(2), nBins);
            end
            xMin = xSpace(startField(i));
            xMax = xSpace(stopField(i));
            posIndX = pos(:, bntConstants.PosX) >= xMin & pos(:, bntConstants.PosX) <= xMax;
            fields(curField).posInd = find(posIndX);
        end

        curField = curField + 1;
    end

    fields(curField:end) = [];
    nFields = length(fields);

    if debug && nFields > 0
        figure, plot(curve);
        v = axis();
        ymin = v(3);
        ymax = v(4);
        color = [1 0 0];
        fieldSize = zeros(1, nFields);
        for i = 1:nFields
            startField = fields(i).col(1);
            stopField = fields(i).col(end);

            p = patch([startField startField stopField stopField], [ymax ymin ymin ymax], color);
            set(p, 'FaceAlpha', 0.5);
            color(1) = color(1) - 0.1;
            color(2) = color(2) + 0.2;
%             color(3) = color(3) + 0.2;
            if color(1) < 0
                color(1) = 0;
            end
            if color(2) > 1
                color(2) = 1;
            end

            fieldSize(i) = fields(i).size;
        end
        xlabel('Bins'), ylabel('Rate');
        titleStr = sprintf('Field(s) with size(s) %s', sprintf('%f, ', fieldSize));
        titleStr(end-1:end) = [];
        title(titleStr);
    end

    % Calculate the spacing using centre of mass (COM). Since centre of mass for a curve is somewhat
    % undefined I'll calculate COM for palte bounded by curve and a function y = 0 for all x.
    % http://tutorial.math.lamar.edu/Classes/CalcII/CenterOfMass.aspx
    distMat = zeros(nFields, 2);

    % store points for debug purposes
    pointsX = zeros(1, nFields);
    pointsY = zeros(1, nFields);

    for i = 1:nFields
        x = fields(i).col(:);
        y = curve(x);
        y = y - min(y); % make it from 0. This gives better estimation

        y1 = x .* y;
        y2 = (y.^2)/2;
        A = trapz(x, y);

        Mx = trapz(x, y1)/A;
        My = min(curve(x)) + trapz(x, y2)/A;

        pointsX(i) = Mx;
        pointsY(i) = My;

        fields(i).x = Mx;
        fields(i).y = My;

        distMat(i, :) = [Mx My];

%             % Raymond's code
%             posMass = sum(mapAxis(x).*curve(x));
%             mass = sum(curve(x));
%             com(i) = posMass/mass;
    end

    if debug
        hold on;
        plot(pointsX, pointsY, 'o');
    end

    if nFields <= 1
        spacing = [];
    else
        D = pdist(distMat);
        spacing = squareform(D);
    end
end