clear all;
close all;
clc
%%
load('ExperimentInformation.mat');
load('GridCellAnalysis.mat');
load('NAT.mat');
%% paremater settings
AngleSmooth=2;
AngleBinsize=3; % cm
SpeedThreadhold=2.5; %cm/mm
MinEventCount=100;
MinSNR=3;
Shuffling=1000;
Shuffling_mininterval=30; %second
%% shuffling
SelectedFrame_raw=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell
SelectedFrame_filtered=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell, filtered by some cratirial
TurningCurve_whole=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
TurningCurve_firsthalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
TurningCurve_secondthalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
TurningCurveStat_whole=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
TurningCurveStat_firsthalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
TurningCurveStat_secondthalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
CellForAnalysis=cell(ExperimentInformation.Session,1);
MVL_shuffled=zeros(ExperimentInformation.TotalCell,Shuffling+3,ExperimentInformation.Session);
Correlation_shuffled=zeros(ExperimentInformation.TotalCell,Shuffling+3,ExperimentInformation.Session);

for j=1:1:ExperimentInformation.Session
    k=1;
    for i=1:1:ExperimentInformation.TotalCell
       SelectedFrame_raw=find(~isnan(NAT{1,j}(:,4*i+10)));
        SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        SelectedFrame_filtered_firstHalf=SelectedFrame_filtered(find(SelectedFrame_filtered<=ExperimentInformation.FrameinEachSession(j)/2));
        SelectedFrame_filtered_secondHalf=SelectedFrame_filtered(find(SelectedFrame_filtered>ExperimentInformation.FrameinEachSession(j)/2));
        HeadDirectionTrain=NAT{1,j}(SelectedFrame_filtered,4);
        HeadDirectionTrain_firstHalf=NAT{1,j}(SelectedFrame_filtered_firstHalf,4);
        HeadDirectionTrain_secondHalf=NAT{1,j}(SelectedFrame_filtered_secondHalf,4);
        Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadholds
        Event_firstHalf=NAT{1,j}(SelectedFrame_filtered_firstHalf,[1 4*i+12]);
        Event_secondHalf=NAT{1,j}(SelectedFrame_filtered_secondHalf,[1 4*i+12]);       
        if length(Event_filtered)>MinEventCount && ExperimentInformation.CellSNR(i)>3 && ~ismember(i,ExperimentInformation.RepeatCell)
            EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
            Event_firstHalf=NAT{1,j}(SelectedFrame_filtered_firstHalf,[1 4*i+12]);
            Event_secondHalf=NAT{1,j}(SelectedFrame_filtered_secondHalf,[1 4*i+12]);
            CellForAnalysis{j,1}(k)=i;
            for n=1:1:Shuffling
                xmin=Shuffling_mininterval*ExperimentInformation.Trackingframerate;
                xmax=size(length(SelectedFrame_filtered),1)-xmin;
                ShiftFrame=round(xmin+rand(1,1)*(xmax-xmin));
                EventTrain_shuffled=circshift(EventTrain,ShiftFrame);
                EventTrain_shuffled(:,1)=EventTrain(:,1);
                EventTrain_shuffled_firstHalf=EventTrain_shuffled(1:size(HeadDirectionTrain_firstHalf,1),:);
                EventTrain_shuffled_secondHalf=EventTrain_shuffled(size(HeadDirectionTrain_firstHalf,1)+1:end,:);
                TurningCurve_shuffled=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain,...
                    EventTrain_shuffled(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                TurnningCurveStastics_shuffled=SpatialTuning_BNT.tcStatistics(TurningCurve_shuffled,AngleBinsize, 50);
                MVL_shuffled(i,n,j)=TurnningCurveStastics_shuffled.r;
                
                TurningCurve_firstHalf_shuffled=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain_firstHalf,...
                    EventTrain_shuffled_firstHalf(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                %                 TurnningCurveStastics_firstHalf=SpatialTuning_BNT.tcStatistics(TurningCurve_firstHalf,AngleBinsize, 50);
                
                TurningCurve_secondHalf_shuffled=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain_secondHalf,...
                    EventTrain_shuffled_secondHalf(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                %                 TurnningCurveStastics_secondHalf=SpatialTuning_BNT.tcStatistics(TurningCurve_secondHalf,AngleBinsize, 50);
                
                Correlation_shuffled(i,n,j)=corr(TurningCurve_firstHalf_shuffled(:,2),TurningCurve_secondHalf_shuffled(:,2));                
                MVL_shuffled(i,Shuffling+1,j)=prctile(MVL_shuffled(i,1:Shuffling,j),95);
                MVL_shuffled(i,Shuffling+2,j)=prctile(MVL_shuffled(i,1:Shuffling,j),99);
                Correlation_shuffled(i,Shuffling+1,j)=prctile(Correlation_shuffled(i,1:Shuffling,j),95);
                Correlation_shuffled(i,Shuffling+2,j)=prctile(Correlation_shuffled(i,1:Shuffling,j),99);               
                TurningCurve_whole{j,i}=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain,...
                    EventTrain(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                TurningCurveStat_whole{j,i}=SpatialTuning_BNT.tcStatistics(TurningCurve_whole{j,i}, AngleBinsize, 50);
                MVL_shuffled(i,Shuffling+3,j)=TurningCurveStat_whole{j,i}.r;            
                TurningCurve_firsthalf{j,i}=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain_firstHalf,...
                    Event_firstHalf(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                TurningCurveStat_firsthalf{j,i}=SpatialTuning_BNT.tcStatistics(TurningCurve_firsthalf{j,i}, AngleBinsize, 50);                
                TurningCurve_secondthalf{j,i}=SpatialTuning_BNT.TurnningCurver_power(HeadDirectionTrain_secondHalf,...
                    Event_secondHalf(:,2),...
                    AngleBinsize,...
                    1/(ExperimentInformation.Trackingframerate/ExperimentInformation.ImagingPlane),...
                    AngleSmooth);
                TurningCurveStat_secondthalf{j,i}=SpatialTuning_BNT.tcStatistics(TurningCurve_secondthalf{j,i},AngleBinsize, 50);                        
                Correlation_shuffled(i,Shuffling+3,j)=corr(TurningCurve_firsthalf{j,i}(:,2),TurningCurve_secondthalf{j,i}(:,2));   
                disp([num2str(i),'-',num2str(n)]);                 
            end
            k=k+1;
        else
        end
    end
end
%% identify HD cells
IsHDCell=cell(ExperimentInformation.Session,1);
for j=1:1:ExperimentInformation.Session
    g=1;
    for i=1:1:ExperimentInformation.TotalCell
            if MVL_shuffled(i,Shuffling+3,j)>MVL_shuffled(i,Shuffling+1,j) && Correlation_shuffled(i,Shuffling+3,j)>Correlation_shuffled(i,Shuffling+1,j)
                IsHDCell{j,1}(g)=i;
                g=g+1;
            else
            end
    end
end
%% save HD cell information
HDCellAnalysis=struct;
HDCellAnalysis.IsHDCell=IsHDCell;
HDCellAnalysis.TurningCurve_whole=TurningCurve_whole;
HDCellAnalysis.TurningCurve_firsthalf=TurningCurve_firsthalf;
HDCellAnalysis.TurningCurve_secondthalf=TurningCurve_secondthalf;
HDCellAnalysis.TurningCurveStat_whole=TurningCurveStat_whole;
HDCellAnalysis.TurningCurveStat_firsthalf=TurningCurveStat_firsthalf;
HDCellAnalysis.TurningCurveStat_secondthalf=TurningCurveStat_secondthalf;
HDCellAnalysis.CellForAnalysis=CellForAnalysis;
HDCellAnalysis.MVL_shuffled=MVL_shuffled;
HDCellAnalysis.Correlation_shuffled=Correlation_shuffled;
HDCellAnalysis.Shuffling=Shuffling;
HDCellAnalysis.AngleSmooth=AngleSmooth;
HDCellAnalysis.AngleBinsize=AngleBinsize; % cm
HDCellAnalysis.SpeedThreadhold=SpeedThreadhold; %cm/mm
HDCellAnalysis.MinEventCount=MinEventCount;
HDCellAnalysis.MinSNR=MinSNR;
HDCellAnalysis.Shuffling_mininterval=Shuffling_mininterval; %second
%% plot  whole-session calcium activity maps, colored by head-direction, for all HD cells
close all
figure
x0=-100;
y0=-100;
width=3440;
height=1200;
set(gcf,'position',[x0,y0,width,height])
j=1;
coluum=20;
S=0;
for k=1:1:size(HDCellAnalysis.IsHDCell{j,1},2)  
    i=HDCellAnalysis.IsHDCell{j,1}(k);
    P=mod(k,coluum);
    if P==0
        P=coluum;
        
    elseif P==1
        S=S+1;
    else
    end
    SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Bestshift=GridCellAnalysis.GridBestShift(i);
    Position=NAT{1,j}(Event_filtered,2:3);
    Position(:,1)=Position(:,1)+ (Bestshift*cos(NAT{1,j}(Event_filtered,4) * pi/180));
    Position(:,2)=Position(:,2)+ (Bestshift*sin(NAT{1,j}(Event_filtered,4) * pi/180));    
    
    Event=NAT{1,j}(Event_filtered,4*i+12);
    Direction=NAT{1,j}(Event_filtered,4);
    Max=max(Event(:));    
    scatter(Position(:,1)+100*(P-1),Position(:,2)-100*(S-1),10*(sqrt(Event./Max)),Direction,'filled','MarkerFaceAlpha',0.8)
    colormap(gca,hsv)
    caxis([0 360]);
%     colorbar;
    ylim([-100*(round(size(HDCellAnalysis.IsHDCell{j,1},2)/coluum)+0.5),50])
    xlim([-50,100*coluum])    
    daspect([1 1 1]);
    hold on    
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);

%% plot  whole-session HD tuning for all HD cells
close all
figure
x0=10;
y0=10;
width=4000;
height=1400;
set(gcf,'position',[x0,y0,width,height])
j=1;
for k=1:1:size(HDCellAnalysis.IsHDCell{j,1},2)  
    i=HDCellAnalysis.IsHDCell{j,1}(k);
    subplot(9,20,k,'align')
    Occupancy=preprocessing.circularSmooth(TurningCurve_whole{j,i}(:,3),AngleSmooth);
    Tuning=TurningCurve_whole{j,i}(:,2);
    MaxTime=max(Occupancy);
    MaxTuning=max(Tuning);
    Occupancy=0.5*Occupancy./MaxTime;
    Tuning=Tuning./MaxTuning;
    H=polarplot([TurningCurve_whole{j,i}(:,1)/180*pi;TurningCurve_whole{j,i}(1,1)/180*pi],[Occupancy;Occupancy(1)],'color',[0.6 0.6 0.6],'LineWidth',1);
    hold on
    H=polarplot([TurningCurve_whole{j,i}(:,1)/180*pi;TurningCurve_whole{j,i}(1,1)/180*pi],[Tuning;Tuning(1)],'color',[0 0 0],'LineWidth',1.5);
    Ax = gca;
    Ax.RTick = [];
    Ax.RTickLabel = [];
    Ax.ThetaTickLabel = [];
    Ax.ThetaTick  = [0 90 180 270];
    Ax.RLim=[0 1];
    AX.FontSize=0;
    Ax.ThetaAxis.Color = [1 1 1];
    AX.RColor = [1 1 1];
    set(gca,'GridColor',[0 0 0]);
    set(gca,'GridAlpha',1);
    set(gca,'RGrid','off');
    set(H,'LineWidth',1);
    AX.RAxis.Visible='off';
%     title(['P-time: ',num2str(MaxTime,'%.2f'),' s',', P-tuning: ',num2str(MaxTuning,'%.2f')])
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);

%%
save ([ExperimentInformation.RawDataAddress,'\HDCellAnalysis.mat'],'HDCellAnalysis','-v7.3');
disp('HD cell analysis was done!');