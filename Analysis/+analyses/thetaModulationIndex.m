% Caclulate theta modulation index
%
% The measure is calculated as described in publication (with small modifications)
% Cacucci et al, Theta-Modulated Place-by-Direction Cellsin the Hippocampal
% Formation in the Rat, Journal of Neuroscience, 2004.
% doi: 10.1523/jneurosci.2635-04.2004
%
% Note that the normalization by trial length is not necessary because it
% doesn't change the resulting value.
%
%  USAGE
%   index = analyses.thetaModulationIndex(spikes)
%   spikes      Vector of spike timestamps in seconds.
%   index       Calculated index, unit-less.
%
function thIndex = thetaModulationIndex(spikes)
    spikes_ms = spikes * 1e3;
    binSize_ms = 5;
    histogramLength_ms = 500;
    sumsMean = zeros([1 floor(histogramLength_ms/binSize_ms)]); % binned count of spikes

    for i = 1:length(spikes_ms)
        curSpike = spikes_ms(i);
        low_edge = curSpike;
        upper_edge = curSpike + histogramLength_ms;
        edges = low_edge:binSize_ms:upper_edge;
        N = histcounts(spikes_ms, edges);
        N(1) = N(1) - 1; % we always have curSpike counted, but we should not count it
        sumsMean = sumsMean + N;
    end
    
    % magic numbers are taken from the publication
    throughBins = [50 70]/binSize_ms; % convert time to bin sizes
    throughBins(1) = floor(throughBins(1));
    throughBins(2) = ceil(throughBins(2));
    peakBins = [100 140]/binSize_ms;
    peakBins(1) = floor(peakBins(1));
    peakBins(2) = ceil(peakBins(2));
    
    through = mean(sumsMean(throughBins(1):throughBins(2)));
    peak = mean(sumsMean(peakBins(1):peakBins(2)));
    % peak - through in order to get mostly positive values
    thIndex = (peak-through) / (through+peak);
end

%% Alternative calculation from Emilio is below. Have not thoroughly compared execution speed/accuracy of the methods.
% function [indx1, altSum] = thetaIndex(ts)

% % h=hist(ts,0:.001:max(ts)+.001);
% fs=1000;
% h=1+round(ts*fs); % convert to ms
% altSum = zeros(length(ts), 500);
% for i = 1:length(ts)
%     altSum(i, :) = arrayfun(@(x) numel(intersect(h(i) + x, h)), 1:500);
% end
% a = arrayfun(@(x) numel(intersect(h+x,h)),1:500);
% A = nanmean(a(100:140));
% B = nanmean(a(50:70));
% indx1 = (A-B) / (A+B);
% histcounts

% sumAsBins = zeros(1, 100);
% altSum1 = sum(altSum);
% for i = 1:100
%     lowEdge = (i-1)*5 + 1;
%     highEdge = lowEdge + 5 - 1;
%     if highEdge > 500
%         highEdge = 500;
%     end
%     sumAsBins(i) = sum(altSum1(lowEdge:highEdge));
% end

% end
