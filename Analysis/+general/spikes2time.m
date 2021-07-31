% Map spike timestamps to some larger time vector
%
% This function return indices of ts that correspond to spikes values.
% The idea is that ts a long signal and spikes are found within it.
% One application is for example to match spikes to LFP signal.
%
%  USAGE
%   spkInd = general.spikes2time(spikes, ts)
%   spikes      Vector of spike timestamps. Units of spikes should be the
%               same as units of ts.
%   ts          Time signal. Vector or matrix. If matrix is provided,
%               the first column is taken as ts.
%   spkInd      Vector of indices in ts that correspond to spike values.
%
function spkInd = spikes2time(spikes, ts)
    if size(ts, 1) > 1 && size(ts, 2) > 1
        ts = ts(:, 1);
    else
        ts = ts(:);
    end
    % make sure it is a column
    spikes = spikes(:);

    spkInd = knnsearch(ts, spikes);
end