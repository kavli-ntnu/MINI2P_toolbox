% Calculate spike mean firing rate
%
% A mean firing rate is obtained by dividing the number of spikes by the duration of a trial/position samples.
% Duration is calculated as the number of tracked position samples multipled by sampling rate. Sampling rate
% is estimated from positions as mean(diff(t)).
%
%  USAGE
%   rate = analyses.meanRate(spikes, positions);
%   spikes      Vector. Spike timestamps.
%   positions   Nx3 or Nx5 ([t x y] or [t x1 y1 x2 y2]) matrix of position samples.
%   rate        Mean firing rate. Units of rate are based on units of time from positions matrix.
%               If time is in seconds, then rate is in Hz.
%
%  EXAMPLE
%
%    data.loadSessions('C:\home\workspace\Li\inputfile-short.cfg');
%    pos = data.getPositions();
%    spikes = data.getSpikeTimes([1 1]);
%    rate = analyses.meanRate(spikes, pos);
%
function rate = meanRate(spikes, positions)
    if isempty(spikes) || isempty(positions)
        rate = 0;
        return;
    end
    
    sampleTime = mode(diff(positions(:, 1)));
    numNans = length(find(isnan(positions(:, 2))));
    duration = (size(positions, 1) - numNans) * sampleTime;
%     duration = (size(positions, 1) ) * sampleTime;
    
    % be sure to match spike to position, so that we do not count extra spikes. Might
    % be usefull for complicated sessions (like trials in a maze or on a linear track).
    if numNans > 0
        spkPos = data.getSpikePositions(spikes, positions);
        spikesLength = size(spkPos, 1);
    else
        spikesLength = length(spikes);
    end
    
%     rate = size(spkPos, 1) / duration;
    rate = spikesLength / duration;
end
