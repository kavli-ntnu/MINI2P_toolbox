% Calculate turning curve
%
% Turning curve resembles a rate map, but calculated based on angular data.
%  USAGE
%
%    tc = analyses.turningCurve(spkhd, poshd, sampleTime, <options>)
%    spkhd          Angles associated with each spike. In degrees. This could be head
%                   directions per spikes, or moving direction.
%    poshd          Vector of angles associated with each position sample (all head directions).
%                   Values should be in degrees.
%    sampleTime     Sample time of recording system, [sec]. If the sampling time is
%                   uniform, then this parameter has no effect on the shape of the curve.
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 1).
%                   Turning curve is created from unsmoothed components, and
%                   only the final result is smoothed.
%     'binWidth'    Bin width in degrees (default = 3).
%    =========================================================================
%    tc             Turning curve matrix of size Nx3. N - number of bins.
%                   1 column: map axis, degrees. Turning curve is based on a histogram, and
%                   map axis consists of bin middle points. For example, if bin
%                   has a range 30-60 degrees, then corresponding map axis value is 45.
%                   2 column: turning curve values.
%                   3 column: trajectory map (how many times head/animal was at some angle). So
%                   it is a histogram binned with the same width as the turning curve.
%
function tc = turningCurve(spkhd, poshd, sampleTime, varargin)
    % Default values
    defaultBinWidth = 3;
    defaultSmooth = 1;
    
    validationVector = @(x) validateattributes(x, {'numeric'}, {'vector', 'nonempty'});
    validationScalar = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonempty'});
    
    inp = inputParser;
    inp.FunctionName = 'analyses.turningCurve';
    addRequired(inp, 'spkhd', validationVector);
    addRequired(inp, 'poshd', validationVector);
    addRequired(inp, 'sampleTime', validationScalar);
    addOptional(inp, 'smooth', defaultSmooth, validationScalar);
    addOptional(inp, 'binWidth', defaultBinWidth, validationScalar);
    
    parse(inp, spkhd, poshd, sampleTime, varargin{:});
    
    binWidth = inp.Results.binWidth;
    smooth = inp.Results.smooth;
    sampleTime = inp.Results.sampleTime;
    
    [spkTurningMap, ~, mapAxis] = general.circHist(spkhd, binWidth);
    allTurningMap = general.circHist(poshd, binWidth);

    
    if ~isequal(size(allTurningMap), size(spkTurningMap))
        spkTurningMap = spkTurningMap';
    end
    
    % Transform from number of sample to amount of time
    allTurningMap = allTurningMap * sampleTime; % assume equal time between samples
    trajectoryMap = allTurningMap;
%     spkTurningMap = general.smooth(spkTurningMap, smooth);
%     allTurningMap = general.smooth(allTurningMap, smooth);
%     trajectoryMap = general.smooth(trajectoryMap, smooth);

    tc(:, 1) = mapAxis;
    tc(:, 2) = spkTurningMap ./ (allTurningMap + eps);
    tc(:, 3) = trajectoryMap;

    tc(:, 2) = general.smooth(tc(:, 2), smooth);
end
