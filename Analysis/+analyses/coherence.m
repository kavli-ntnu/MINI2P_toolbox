% Calculate spatial coherence of a rate map
%
% The value is calculated according to RU Muller, JL Kubie "The firing of hippocampal place
% cells predicts the future position of freely moving rats", Journal of Neuroscience, 1 December 1989,
% 9(12):4101-4110. The paper doesn't provide information about how to deal with border values
% which do not have 8 well-defined neighbours. This function uses zero-padding technique.
%
%  USAGE
%   z = analyses.coherence(map, <options>)
%
%   map            2D rate map, could have NaNs. NaNs are replaced with 0 for the calculation.
%   <options>      optional list of property-value pairs (see table below)
%
%   =============================================================================
%    Properties    Values
%   -----------------------------------------------------------------------------
%    'normalize'    normalize the result (default = 'on'). Other value is 'off'.
%   =============================================================================
%
%   z               Coherence value
%
function z = coherence(map, varargin)
    inp = inputParser;
    defaultNormalization = 'on';

    checkNormaliztion = @(x) strcmpi(x, 'on') || strcmpi(x, 'off');
    checkMap = @(x) ismatrix(x) && ~isempty(x);

    addRequired(inp, 'map', checkMap);
    addParamValue(inp, 'normalize', defaultNormalization, checkNormaliztion);

    parse(inp, map, varargin{:});

    % get parsed arguments
    doNormalization = strcmpi(inp.Results.normalize, 'on');

    % averaging kernel, 1/8
    K = [0.125 0.125 0.125;
         0.125     0 0.125;
         0.125 0.125 0.125];

    nanIdx = isnan(map);
    if any(nanIdx(:))
        map(nanIdx) = 0;
    end

    avgMap = conv2(map, K, 'same');
    avgMapLinear = reshape(avgMap', numel(avgMap), 1); % reshape linearly row-based
    mapLinear = reshape(map', numel(map), 1);

    z = corr(mapLinear, avgMapLinear);
    if doNormalization
        z = atanh(z);
    end
end