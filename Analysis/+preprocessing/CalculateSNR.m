function SNR=CalculateSNR(RawTrancient,signalTime,restingTime)

RestingPeriod=RawTrancient(restingTime);
SignalPeriod=RawTrancient(signalTime);

Peaks=SignalPeriod(find(SignalPeriod>=prctile(SignalPeriod,90)));  % definition of peaks: large than 99th intensity, large than 2x std of the traces
% Baseline=RawTrancient(find(RawTrancient<prctile(SignalPeriod,1)));  % 
Baseline=RestingPeriod;
Signal=median(Peaks)-median(Baseline);
% FrameDifference=abs(RestingPeriod-circshift(RestingPeriod,1));
% noise=median(FrameDifference); % the extraction of the average noise level is from the paper "A deep learning toolbox for noise-optimized, generalized spike inference from calcium imaging data",https://www.biorxiv.org/content/10.1101/2020.08.31.272450v1
noise=std(RawTrancient); % t
disp(['Noise level is ', num2str(noise),'; Signal level is ',num2str(Signal),]);
SNR=Signal/noise;
if isnan(SNR)
    SNR=0;
else
end
end


