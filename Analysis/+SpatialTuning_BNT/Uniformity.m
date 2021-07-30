function EOC_Uniformity=Uniformity(Spike,Distance,Angle,DistanceRange,degSamp,framerate)
% i=22;
% j=2;
% SpikeTrain=FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,j}(:,4*i+12);
% SelectFrame=find(~isnan(SpikeTrain).*(FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,j}(:,5)>FullNeuronBehaviorDataSet.SpeedThreadhold)==1);
% Distance=FullNeuronBehaviorDataSet.RToObject{1,j};
% Distance=Distance(SelectFrame);
% Angle=FullNeuronBehaviorDataSet.AngleToObject{1,j};  
% Angle=Angle(SelectFrame);
% DistanceRange=FullNeuronBehaviorDataSet.EOC_TurningStastics{j,i}.Distance_FWHM;  
% Spike=SpikeTrain(SelectFrame);
% degSamp=30;
% framerate=FullNeuronBehaviorDataSet.FrameRate./double(FullNeuronBehaviorDataSet.ImagingPlane);
NumAngularBins=round(360/degSamp)+1; 
AngularBins=linspace(0,360,NumAngularBins);
EOC_Uniformity=struct;
EOC_Uniformity.Occ_infield=zeros(NumAngularBins,1);
EOC_Uniformity.Occ_outfield=zeros(NumAngularBins,1);
EOC_Uniformity.Spike_infield=zeros(NumAngularBins,1);
EOC_Uniformity.Spike_outfield=zeros(NumAngularBins,1);
for i=1:1:size(Distance,1)
    if Distance(i)>=DistanceRange(1)&& Distance(i)<=DistanceRange(2)
        for j=1:1:NumAngularBins-1
            if Angle(i)>=AngularBins(j)&& Angle(i)<AngularBins(j+1)
              EOC_Uniformity.Occ_infield(j)=EOC_Uniformity.Occ_infield(j)+1;
              EOC_Uniformity.Spike_infield(j)=EOC_Uniformity.Spike_infield(j)+Spike(i);
            else
            end
        end
    else
        for j=1:1:NumAngularBins-1
            if Angle(i)>=AngularBins(j)&& Angle(i)<AngularBins(j+1)
              EOC_Uniformity.Occ_outfield(j)=EOC_Uniformity.Occ_outfield(j)+1;
              EOC_Uniformity.Spike_outfield(j)=EOC_Uniformity.Spike_outfield(j)+Spike(i);
            else
            end
        end
    end
end
sampleTime=1/framerate;
EOC_Uniformity.Occ_outfield=EOC_Uniformity.Occ_outfield.*sampleTime;   
EOC_Uniformity.Occ_outfield=EOC_Uniformity.Occ_outfield(1:end-1);
EOC_Uniformity.Occ_infield=EOC_Uniformity.Occ_infield.*sampleTime;   
EOC_Uniformity.Occ_infield=EOC_Uniformity.Occ_infield(1:end-1);
EOC_Uniformity.Spike_outfield=EOC_Uniformity.Spike_outfield(1:end-1);
EOC_Uniformity.Spike_infield=EOC_Uniformity.Spike_infield(1:end-1);
EOC_Uniformity.Rate_infield=EOC_Uniformity.Spike_infield./EOC_Uniformity.Occ_infield;
EOC_Uniformity.Rate_outfield=EOC_Uniformity.Spike_outfield./EOC_Uniformity.Occ_outfield;
EOC_Uniformity.Uniformity=(EOC_Uniformity.Rate_infield-EOC_Uniformity.Rate_outfield)>0;
EOC_Uniformity.UniformityScore=sum(EOC_Uniformity.Uniformity)./size(EOC_Uniformity.Uniformity,1);


end

