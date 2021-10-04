% Calculate speed score(s)
%
% Speed score is a correlation between cell firing rate and animal's speed.
% This is adaptation of code from Emilio Kropff, ekropff@leloir.org.ar
% The script calculates two types of speed score refered to as 'score15' and
% 'score16'. Score15 is a speed score as appeared in paper doi:10.1038/nature14622.
% Note that in order to match Emilio's original code, you should only use a low
% boundary of speed filter.
% Score16 is a new version of the speed score from Emilio.
% The difference between these score lies in the way how areas with invalid speed
% are processed.
% Score15 replaces samples of firing rate that fall to invalid speed with NaNs
% and then smooths the firing rate.
% Score16 smooths the instanteneous firing rate as is, then takes only valid
% speed indices into account during the correlation step.
%
%  USAGE
%   speedScores = analyses.speedScore(speed, firingRate, span, <options>)
%   speed       Nx1 vector of speed measurements.
%   firingRate  Nx1 vector of instanteneous firing rate, unsmoothed!
%   span        Gaussian filter span/width in samples. Gaussian filter is used
%               to smooth firing rate. Original paper used 0.4 sec.
%   <options>   Optional list of property-value pairs (see table below)
%
%   ============================================================================
%    Properties           Values
%   ----------------------------------------------------------------------------
%    speedIndices          Nx1 logical vector that defines valid speed samples.
%                          This can be result of applying speed filter to speed.
%                          Default value: validSpeedIndices = true(size(speed));
%   ============================================================================
%   speedScores  Speed score values [score15 score16].
%
%  SEE
%
%    See also general.smoothGauss.
%
function scores = speedScore(speed, firingRate, span, varargin)
    inp = inputParser;

    defaultSpeedIndices = true(size(speed));

    checkSpeedIndices = @(x) isequal(size(x), size(speed));

    addRequired(inp, 'speed');
    addRequired(inp, 'firingRate');
    addRequired(inp, 'span');
    addParameter(inp, 'speedIndices', defaultSpeedIndices, checkSpeedIndices);

    if isempty(speed) || isempty(firingRate) || span <= 1
        scores = [nan nan];
        return;
    end

    parse(inp, speed, firingRate, span, varargin{:});
    goodSpeed = inp.Results.speedIndices;
    badSpeed = ~goodSpeed;

    % score15
    frWithNan = firingRate;
    frWithNan(badSpeed) = NaN;
    frWithNanSmoothed = general.smoothGauss(frWithNan, span);
    scores(1) = corr(speed(goodSpeed), frWithNanSmoothed(goodSpeed), 'rows', 'pairwise','type','Spearman');

    % score16
    firingRateSmoothed = general.smoothGauss(firingRate, span);
    scores(2) = corr(speed(goodSpeed), firingRateSmoothed(goodSpeed), 'rows', 'pairwise','type','Spearman');
end
