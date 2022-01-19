% Detect spikes in an EEG signal
%
% Detect spike by means of a threshold.
%
%  USAGE
%   spikes = analyses.detectSpikes(eeg, sampleFreq, spikeWidth_sec, <options>)
%   eeg             Data vector with spikes
%   sampleFreq      Data sampling frequency in Hz
%   spikeWidth_sec  Width of a single spike in seconds (default 2e-3)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bandpass'    [fl fh], bandpass frequencies. Values in Hz that is used to
%                   bandpass eeg signal (default is [600 3000]).
%     'threshold'   Threshold value. Either a single value [th] or vector of two
%                   elements [nt st].
%                   A single value th means that the spike threshold is calculated as:
%                       spikeTh = th;
%                   th most likely is an absolute value in mV. For example th = 70.
%                   In case of [nt st] values provided, the threshold for spike detection
%                   is calculated as:
%                       noiseTh = nt * std(eeg);
%                       spikeTh = st * noiseTh;
%                   Default value is [2 1.5].
%     'polarity'    {'negative', 'positive', 'both'}, string that defines what type of spike
%                   events to detect. 'Negative' means that only negative spikes are detected,
%                   'positive' means that only positive spikes are used, 'both' stays for
%                   both types of spikes. Negative and positive differ in the direction of spike
%                   rise. Positive being towards the higher values and negative towards smaller
%                   values. Default value is 'both'.
%    =========================================================================
%   spikes          Vector of spike timestamps in sec.
%   spikeWaves      NxM matrix with spike waveforms. N is the number of detected spikes.
%                   M depends on the sampling frequency. 1 millisecond of data is extracted
%                   for each spike (200 µs before and 800 µs after the threshold event).
%                   M/sampleFreq = 1 ms.
%
function [spikes, spikeWaves] = detectSpikes(eeg, sampleFreq, spikeWidth_sec, varargin)
    time = [1:length(eeg)] / sampleFreq;
    spikeWidth = spikeWidth_sec * sampleFreq; % spike width in number of points

    spikeWaves = [];

    inp = inputParser;
    defaultBandpass = [600 3000];
    defaultThreshold = [2 1.5];
    defaultPolarity = 'both';

    checkBandpass = @(x) length(x) == 2;
    checkThreshold = @(x) helpers.isdvector(x, '>0');
    checkPolarity = @(x) helpers.isstring(lower(x), 'negative', 'positive', 'both');

    addRequired(inp, 'eeg');
    addRequired(inp, 'sampleFreq');
    addRequired(inp, 'spikeWidth_sec');
    addParamValue(inp, 'bandpass', defaultBandpass, checkBandpass);
    addParamValue(inp, 'threshold', defaultThreshold, checkThreshold);
    addParamValue(inp, 'polarity', defaultPolarity, checkPolarity);

    parse(inp, eeg, sampleFreq, spikeWidth_sec, varargin{:});

    bandpass = inp.Results.bandpass;
    thValue = inp.Results.threshold;
    doStd = length(thValue) == 2;
    polarity = lower(inp.Results.polarity);

    data2 = general.fastFftBandpass(eeg, sampleFreq, bandpass(1), bandpass(2));

    % thresholds
    if doStd
        noiseTh = thValue(1) * std(data2);
        spikeTh = thValue(2) * noiseTh;
    else
        spikeTh = thValue(1);
    end

    % threshold at a set level
    switch polarity
        case 'negative'
            possibleSpikes = find(data2 < -spikeTh);
            
        case 'positive'
            possibleSpikes = find(data2 > spikeTh);
            
        case 'both'
            possibleSpikes = find(data2 > spikeTh | data2 < -spikeTh);
    end

    % clean possible spikes to get rid of duplicates

    % find one spike in each group
%     spikeWidthMax = 40e-3 * sampleFreq;
    diffSpikes = diff(possibleSpikes);
%     dumspikes = possibleSpikes(diffSpikes > spikeWidth & diffSpikes < spikeWidthMax);
    dumspikes = possibleSpikes(diffSpikes > spikeWidth);

    % use that one spike as an index to find
    %   the maximum spike in the group
    spikes = zeros(1, length(dumspikes)-1);
    for i=1:length(dumspikes)-1
        maxingroup=find(max( abs( data2(dumspikes(i)+1:dumspikes(i+1)) ))...
            == abs(data2(dumspikes(i)+1:dumspikes(i+1))));
        spikes(i) = time(maxingroup+dumspikes(i));
    end

    if nargout > 1
        % 200 micro-sec before and 800 micro-sec after
        beginOffset_sec = 200e-6;
        endOffset_sec = 800e-6;
        waveBegins_sec = spikes - beginOffset_sec;
        waveEnds_sec = spikes + endOffset_sec;
        waveBegins = round(waveBegins_sec * sampleFreq);
        waveEnds = round(waveEnds_sec * sampleFreq);
        waveWidth = waveEnds(1) - waveBegins(1) + 1;

        spikeWaves = nan(length(spikes), waveWidth);
        for i = 1:length(spikes)
            spikeWaves(i, :) = eeg(waveBegins(i):waveEnds(i));
        end
    end
end