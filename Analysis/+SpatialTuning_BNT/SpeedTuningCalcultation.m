function SpeedTuning=SpeedTuningCalcultation(SpeedTrain,SpikeTrain,FrameRate,SpeedRange,MinSpanTime,extend)
% FrameRate=FullNeuronBehaviorDataSet.FrameRate./double(FullNeuronBehaviorDataSet.ImagingPlane);
% Frameselect=find(~isnan(FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,2}(:,12+4*i)));
% SpeedTrain=FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,2}(Frameselect,5);
% SpikeTrain=FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,2}(Frameselect,4*i+12);
% SpeedRange=[0 4 8 12];
% MinSpanTime=10;

SpeedTuning=struct;
SpeedTuning.Rate=zeros(size(SpeedRange));
SpeedTuning.Occ=zeros(size(SpeedRange));
SpeedTuning.SpikeCount=zeros(size(SpeedRange));
SpeedTuning.SpeedSelectivity=0;
SpeedTuning.SpeedScore=0;
SpeedTuning.PeakRate=0;
SpeedTuning.PeakSpeed=1;
SpeedTuning.SpeedRange=SpeedRange;
SpeedRangeExpend=[SpeedRange,extend];
for j=1:1:length(SpeedTrain)
    for k=1:1:size(SpeedRangeExpend,2)-1
        if SpeedTrain(j)>=SpeedRangeExpend(k)&&SpeedTrain(j)<SpeedRangeExpend(k+1)
            SpeedTuning.Occ(k)=SpeedTuning.Occ(k)+1;
            SpeedTuning.SpikeCount(k)=SpeedTuning.SpikeCount(k)+SpikeTrain(j);
        else
        end
    end
end
SpeedTuning.Occ=SpeedTuning.Occ./FrameRate;% frame to second;
SpeedTuning.Occ(SpeedTuning.Occ==0)=0.0000001;
SpeedTuning.Rate=SpeedTuning.SpikeCount./SpeedTuning.Occ;
% SpeedTuning.Rate=smooth(SpeedTuning.Rate,speedsmooth);
SpeedTuning.Rate(SpeedTuning.Occ<MinSpanTime)=nan;

[SpeedTuning.PeakRate,SpeedTuning.PeakSpeed]=max(SpeedTuning.Rate(~isnan(SpeedTuning.Rate)));
% SpeedTuning.PeakSpeed=SpeedRange(SpeedTuning.PeakSpeed);
MinRate=min(SpeedTuning.Rate(~isnan(SpeedTuning.Rate)));
SpeedTuning.SpeedSelectivity=(SpeedTuning.PeakRate-MinRate)./(SpeedTuning.PeakRate+MinRate);

% Con1=find((SpeedRange<=ConditionSeperation).*(~isnan(SpeedTuning.Rate))==1);
% Con2=find((SpeedRange>ConditionSeperation).*(~isnan(SpeedTuning.Rate))==1);
% Condtion1=mean(SpeedTuning.Rate(Con1));
% Condtion2=mean(SpeedTuning.Rate(Con2));
% SpeedTuning.SpeedScore=(Condtion2-Condtion1)./(Condtion1+Condtion2);

end








% end