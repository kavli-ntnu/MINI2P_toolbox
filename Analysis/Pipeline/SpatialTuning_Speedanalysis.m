close all;
%%
load('ExperimentInformation.mat');
load('NAT.mat');
%% paremater settings
SpeedBinning=1.5; %cm/s
MinSpeed=2.5;
MaxSpeed=16;
SpeedRange=[MinSpeed:SpeedBinning:MaxSpeed]; %cm0
MinSpanTime=10;%second
ShuffleNum=200;
Shuffling_mininterval=30;
SpeedSmooth=round(0.25*ExperimentInformation.Trackingframerate);
MinEventCount=100;

%%
SpeedTunning=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); % tuning curve
Speed_KStest=zeros(ExperimentInformation.Session,ExperimentInformation.TotalCell,ShuffleNum+3); % K-S value for shuffles (1 to ShuffleNum, 95th, 5th, and observed data)
Speed_Spearman=zeros(ExperimentInformation.Session,ExperimentInformation.TotalCell,ShuffleNum+3); % Spearman for shuffles (1 to ShuffleNum, 95th, 5th, and observed data)

%% K-S analysis and Spearmean test
% generate speed tuning curve with finner binning

for j=1:1:ExperimentInformation.Session
    k=1;
    for i=1:1:ExperimentInformation.TotalCell
        
        Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>MinSpeed));% filter out the frames with speed threadholds
       if length(Event_filtered)>MinEventCount && ExperimentInformation.CellSNR(i)>3 && ~ismember(i,ExperimentInformation.RepeatCell)
       SelectedFrame_raw=find(~isnan(NAT{1,j}(:,12+4*i))); %all the position are select
       SpeedTrain=NAT{1,1}(SelectedFrame_raw,5);
       SpeedTrain_smooth=general.smoothGauss(SpeedTrain,SpeedSmooth); 
       SelectedFrame_filtered=logical((SpeedTrain_smooth>=MinSpeed).*(SpeedTrain_smooth<=MaxSpeed+SpeedBinning));
       SpeedTrain_filtered=SpeedTrain_smooth(SelectedFrame_filtered);
       EventTrain_original=NAT{1,j}(SelectedFrame_raw,4*i+12);
       EventTrain_smooth=general.smoothGauss(EventTrain_original,SpeedSmooth);
       EventTrain_filtered=EventTrain_smooth(SelectedFrame_filtered);    
       SpeedTunning{j,i}=SpatialTuning_BNT.SpeedTuningCalcultation(SpeedTrain_filtered,EventTrain_filtered,...
                         ExperimentInformation.Trackingframerate./double(ExperimentInformation.ImagingPlane),...
                         SpeedRange,MinSpanTime,...
                         MaxSpeed+SpeedBinning);
       SpeedTunning{j,i}.Rate=smoothdata(SpeedTunning{k,1}.Rate,'lowess',3);
       SpeedCount=round(SpeedTunning{j,i}.Rate./max(SpeedTunning{j,i}.Rate)*100);
       KS_SpeedBinVector{k,1}=zeros(sum(SpeedCount(~isnan(SpeedCount))),1);
       m=1;
              for p=1:1:length(SpeedCount)
                  if ~isnan(SpeedCount(p))
                    for q=1:1:SpeedCount(p)
                            KS_SpeedBinVector{k,1}(m,1)=SpeedTunning{k,1}.SpeedRange(p);
                            m=m+1;
                    end
                  else
                  end
              end   
      for n=1:1:ShuffleNum              
          xmin=Shuffling_mininterval*ExperimentInformation.Trackingframerate;
          xmax=size(SelectedFrame_filtered,1)-xmin;
          ShiftFrame=round(xmin+rand(1,1)*(xmax-xmin));
          SpikeTrain_raw_shuffled=circshift(EventTrain_filtered,ShiftFrame);           
          Speedtuning_shuffled=SpatialTuning_BNT.SpeedTuningCalcultation(SpeedTrain_filtered,SpikeTrain_raw_shuffled,...
                                ExperimentInformation.Trackingframerate./double(ExperimentInformation.ImagingPlane),...
                                SpeedTunning{k,1}.SpeedRange,MinSpanTime,17.5);  
          RateSmoothed=smoothdata(Speedtuning_shuffled.Rate,'lowess',3);
          SpeedCount=round(RateSmoothed./max(SpeedTunning{k,1}.Rate)*100);         

          % calculate the KS value for shuffles
          KS_SpeedBinVector_shuffle=zeros(sum(SpeedCount(~isnan(SpeedCount))),1);
          m=1;
          for p=1:1:length(SpeedCount)
              if ~isnan(SpeedCount(p))
                  for q=1:1:SpeedCount(p)
                      KS_SpeedBinVector_shuffle(m,1)=Speedtuning_shuffled.SpeedRange(p);
                      m=m+1;
                  end
              else
              end
          end
          [~,~,Speed_KStest(j,i,n)]=kstest2(KS_SpeedBinVector_shuffle,SpeedTunning{k,1}.SpeedRange);
          % calculate the spearman correlation
          speedScore_tem=SpatialTuning_BNT.speedScore(SpeedTrain_filtered,SpikeTrain_raw_shuffled(:,2),...
                         ExperimentInformation.Trackingframerate./double(ExperimentInformation.ImagingPlane));
          Speed_Spearman(j,i,n)=speedScore_tem(1);  
     
      end
      Speed_KStest(j,i,ShuffleNum+1)=prctile(Speed_KStest(j,i,1:ShuffleNum),95);
      Speed_KStest(j,i,ShuffleNum+2)=prctile(Speed_KStest(j,i,1:ShuffleNum),5);
      [~,~,Speed_KStest(k,ShuffleNum+3)]=kstest2(KS_SpeedBinVector{k,1},SpeedTunning{k,1}.SpeedRange);   
      
      Speed_Spearman(j,i,ShuffleNum+1)=prctile(Speed_Spearman(j,i,1:ShuffleNum),95);
      Speed_Spearman(j,i,ShuffleNum+2)=prctile(Speed_Spearman(j,i,1:ShuffleNum),5);      
      speedScore_tem=SpatialTuning_BNT.speedScore(SpeedTrain_filtered,EventTrain_filtered(:,2),...
                                                  ExperimentInformation.Trackingframerate./double(ExperimentInformation.ImagingPlane));
      Speed_Spearman(j,i,ShuffleNum+3)=speedScore_tem(1);
       else
       end 
    end
end

%% identify speed cells
IsPosSpeedCell=cell(ExperimentInformation.Session,1);
IsNagSpeedCell=cell(ExperimentInformation.Session,1);
IsNonlinearSpeedCell=cell(ExperimentInformation.Session,1);
for j=1:1:ExperimentInformation.Session
    g=1;
    k=1;
    h=1;
    for i=1:1:ExperimentInformation.TotalCell
            if Speed_KStest(j,i,Shuffling+3,j)>Speed_KStest(j,i,Shuffling+1,j)
                if Speed_Spearman(j,i,ShuffleNum+3)>Speed_Spearman(j,i,ShuffleNum+1)
                    IsPosSpeedCell{j,1}(g)=i;
                    g=g+1;
                elseif Speed_Spearman(j,i,ShuffleNum+3)<Speed_Spearman(j,i,ShuffleNum+2)
                    IsNagSpeedCell{j,1}(k)=i;
                    k=k+1;
                else
                    IsNonlinearSpeedCell{j,1}(h)=i;
                    h=h+1;       
                end
            else    
            end
    end
end    
    