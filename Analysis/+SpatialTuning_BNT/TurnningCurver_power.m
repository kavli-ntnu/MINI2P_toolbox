function tc=TurnningCurver_power(angle,amplitute,binwidth,sampleTime,smoothfactor)

% input:
% amplitute: the calcium signal in each frame;
% angle: the angle in each frame;
% the binwidth to make the histogrom;

% amplitute must have the same length as angle 

angle = mod(angle, 360); % Make sure 0 <= theta <= 2*pi
numBins = ceil(360 / binwidth);
x = (0:numBins-1) * 360/numBins + 180./numBins;
tc(:, 1)= x;
%     edges = sort(mod([(x(2:end) + x(1:end-1))/2 (x(end) + x(1) + twoPi)/2], twoPi));
%     edges = [edges edges(1) + twoPi];

[timeperbin, edge] = histcounts(angle(:,1), numBins(:,1), 'binlimits', [0 360]);
Amplituteperbin=zeros(numBins,1);

for i=1:1:size(angle,1)
    if ceil(angle(i)/binwidth)==0
        Amplituteperbin(1)=Amplituteperbin(1)+amplitute(i);
    else
        Amplituteperbin(ceil(angle(i)/binwidth))=Amplituteperbin(ceil(angle(i)/binwidth))+amplitute(i);
    end
end
timeperbin=timeperbin*sampleTime;
tc(:, 2)=Amplituteperbin./(timeperbin'+eps);
tc(:, 3)=timeperbin';
tc(:, 2)=general.smooth(tc(:, 2), smoothfactor);

end