% Remove unwanted periods in positions, recreate time and move spike times accordingly
%
% This function can be used to remove some of position samples keeping spike timestamps
% in sync. One application is to remove positions that correspond to slow speed.
%
%  USAGE
%   [newPos, spkInd, newCounter] = general.trimPositionAndSpikes(pos, spikeCounter, <options>)
%   pos             NxM matrix of position samples in form [t x y ...]. The very first
%                   column must contain time.
%   spikeCounter    Nx1 vector, i.e. it contains the same number of samples as pos.
%                   Each element in spikeCounter indicates how many spikes fall to that
%                   particular position sample. spikeCounter(3) = 2 indicates that 2 spikes
%                   correspond to position sample 3. And spikeCounter(5) = 0 indicates that
%                   no spikes correspond to position sample 5.
%   <options>       Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties     Values
%   ----------------------------------------------------------------------------------------------
%    'trim'         Nx1 vector. This is a logical vector that indicates which position samples should
%                   be removed. If not provided, then trim = isnan(pos(:, 2)), i.e. trim vector
%                   will be indices of NaN values in second column of pos matrix.
%   ==============================================================================================
%   newPos          JxM matrix of trimmed positions in form [t x y ...]. Time column is recreated
%                   to be linear with sampling rate extracted from pos.
%   spkInd          Linear indices of spikes in newPos matrix, size Kx1.
%   newCounter      Trimmed version of spikeCounter, size Jx1.
%
function [newPos, spkInd, spikeCounter] = trimPositionAndSpikes(pos, spikeCounter, varargin)
    inp = inputParser;
    deafaultTrim = [];

    addRequired(inp, 'pos');
    addRequired(inp, 'spikeCounter');
    addParameter(inp, 'trim', deafaultTrim);
    parse(inp, pos, spikeCounter, varargin{:});
    badPos = inp.Results.trim;
    if isempty(badPos)
        badPos = isnan(pos(:, 2));
    else
        if numel(badPos) ~= size(pos, 1)
            error('BNT:arg:length', 'Parameter ''trim'' should have the same number of elements as number of rows in ''pos'' matrix');
        end
    end

    sampleTime = mean(diff(pos(:, bntConstants.PosT)));
    if isnan(sampleTime)
        sampleTime = data.sampleTime('sec');
    end

    newPos = pos;
    newPos(badPos, :) = [];
    spikeCounter(badPos, :) = [];

    % make new timestamps
    newPos(:, 1) = 0:sampleTime:(size(newPos, 1) - 1) * sampleTime;
    %spikes = repelem(newPos(:, bntConstants.PosT), spikeCounter);
    spkInd = repelem(1:size(newPos, 1), spikeCounter);
end