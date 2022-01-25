% Calculate spike instantaneous firing rate
%
% Instantaneous firing rate is calculated per each position sample and is expressed
% in Hz.
%
%  USAGE
%       firingRate = analyses.instantRate(spikes, pos)
%       spikes      Nx1 vector of spike timestamps. Timestamps should be already matched
%                   to positions! See data.getSpikePositions.
%       pos         Mx3 or Mx5 matrix of position samples.
%       firingRate  Mx1 vector with firing rate expressed in Hz. Or more generally,
%                   units of firing rate depend on units of position timestamps.
%                   If positions are in seconds, then firing rate is in Hz.
%
function firingRate = instantRate(spikes, pos)
    n = histc(spikes, pos(:, 1)); % count in number of occurrences
    dt = diff(pos(:, 1));
    if isempty(dt)
        dt = 0;
    end
    dt(end+1) = dt(end);

    firingRate = [];
    try
        firingRate = n ./ dt;
    catch
        if ~isequal(size(n), size(dt)) && isequal(size(n'), size(dt))
            firingRate = n' ./ dt;
        end
    end
end