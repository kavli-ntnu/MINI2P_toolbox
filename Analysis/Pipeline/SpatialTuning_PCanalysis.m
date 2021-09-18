close all;
% paremater settings
MapSmooth=3;
MapBinsize=2.5; % cm
SpeedThreadhold=2.5; %cm/mm
MinEventCount=100;
MinSNR=3;
MinBins_field=9;
MinPeek_field=0.01;
Cutoff_field=0.5;
Shuffling=1000;
Shuffling_mininterval=30;
MinTime=0.1;
MaxFieldsize=625; %16x16
%% shuffling
SelectedFrame_raw=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell
SelectedFrame_filtered=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell, filtered by some cratirial
ActivityMap=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
ActivityMap_firstHalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
ActivityMap_secondHalf=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
ActivityMapStats=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %%filtered the map with 1) speed threadhold, moving condition or other critirial.
ActivityMapCrossCorrelation=zeros(1,ExperimentInformation.TotalCell);
CellForAnalysis=cell(ExperimentInformation.Session,1);
SI_shuffled=zeros(ExperimentInformation.TotalCell,Shuffling+3,ExperimentInformation.Session);
MapCorrelation_shuffled=zeros(ExperimentInformation.TotalCell,Shuffling+3,ExperimentInformation.Session);
for j=1:1:ExperimentInformation.Session
    k=1;
    for i=1:1:ExperimentInformation.TotalCell
        SelectedFrame_raw=find(~isnan(NAT{1,j}(:,4*i+10)));
        SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        
        if length(Event_filtered)>MinEventCount && ExperimentInformation.CellSNR(i)>3 && ~ismember(i,ExperimentInformation.RepeatCell)
            SelectedFrame_firstHalf=SelectedFrame_filtered(find(SelectedFrame_filtered<=ExperimentInformation.FrameinEachSession/2));
            SelectedFrame_secondHalf=SelectedFrame_filtered(find(SelectedFrame_filtered>ExperimentInformation.FrameinEachSession/2));
            PositionTrain=NAT{1,j}(SelectedFrame_filtered,1:3);
            PositionTrain_firstHalf=NAT{1,j}(SelectedFrame_firstHalf,1:3);
            PositionTrain_secondHalf=NAT{1,j}(SelectedFrame_secondHalf,1:3);
            EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
            EventTrain_firstHalf=NAT{1,j}(SelectedFrame_firstHalf,[1 4*i+12]);
            EventTrain_secondHalf=NAT{1,j}(SelectedFrame_secondHalf,[1 4*i+12]);
            CellForAnalysis{j,1}(k)=i;
            ActivityMap{j,i}=analyses.map(PositionTrain,EventTrain,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',[-32 32 -32 32]);
            ActivityMap_firstHalf{j,i}=analyses.map(PositionTrain_firstHalf,EventTrain_firstHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',0.1,'limits',[-32 32 -32 32]);
            ActivityMap_secondHalf{j,i}=analyses.map(PositionTrain_secondHalf,EventTrain_secondHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',0.1,'limits',[-32 32 -32 32]);
            ActivityMapStats{j,i}=analyses.mapStatsPDF(ActivityMap{j,i});
            ActivityMapCrossCorrelation(j,i)=analyses.spatialCrossCorrelation(ActivityMap_firstHalf{j,i}.z, ActivityMap_secondHalf{j,i}.z);
            
            for n=1:1:Shuffling
                xmin=Shuffling_mininterval*ExperimentInformation.FrameRate;
                xmax=size(length(SelectedFrame_filtered),1)-xmin;
                ShiftFrame=round(xmin+rand(1,1)*(xmax-xmin));
                EventTrain_shuffled=circshift(EventTrain,ShiftFrame);
                EventTrain_shuffled_firstHalf=EventTrain_shuffled(1:size(PositionTrain_firstHalf,1),:);
                EventTrain_shuffled_secondHalf=EventTrain_shuffled(size(PositionTrain_firstHalf,1)+1:end,:);
                ActivityMap_shuffled=analyses.map(PositionTrain,EventTrain_shuffled,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',[-32 32 -32 32]);
                ActivityMapStats_shuffled=analyses.mapStatsPDF(ActivityMap_shuffled);
                SI_shuffled(i,n,j)=ActivityMapStats_shuffled.content;
                Map_shuffled_firstHalf=analyses.map(PositionTrain_firstHalf,EventTrain_shuffled_firstHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime/2,'limits',[-32 32 -32 32]);
                Map_shuffled_secondHalf=analyses.map(PositionTrain_secondHalf,EventTrain_shuffled_secondHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime/2,'limits',[-32 32 -32 32]);
                MapCorrelation_shuffled(i,n,j)=analyses.spatialCrossCorrelation(Map_shuffled_firstHalf.z, Map_shuffled_secondHalf.z);
            end
            k=k+1;
        else
        end
        i
    end
    
end
%% identify PC candidate
IsPCCell_candidate=cell(ExperimentInformation.Session,1);

for j=1:1:ExperimentInformation.Session
    for i=1:1:ExperimentInformation.TotalCell
        if ~isempty(ActivityMapStats{j,i})
            SI_shuffled(i,Shuffling+1,j)=prctile(SI_shuffled(i,1:Shuffling,j),95);
            SI_shuffled(i,Shuffling+2,j)=prctile(SI_shuffled(i,1:Shuffling,j),99);
            SI_shuffled(i,Shuffling+3,j)=ActivityMapStats{j,i}.content;
            MapCorrelation_shuffled(i,Shuffling+1,j)=prctile(MapCorrelation_shuffled(i,1:Shuffling,j),95);
            MapCorrelation_shuffled(i,Shuffling+2,j)=prctile(MapCorrelation_shuffled(i,1:Shuffling,j),99);
            MapCorrelation_shuffled(i,Shuffling+3,j)=ActivityMapCrossCorrelation(j,i);
        else
        end
    end
end
for j=1:1:ExperimentInformation.Session
    g=1;
    for k=1:1:ExperimentInformation.TotalCell
        if SI_shuffled(k,Shuffling+3,j)>SI_shuffled(k,Shuffling+1,j)
            if MapCorrelation_shuffled(k,Shuffling+3,j)>MapCorrelation_shuffled(k,Shuffling+1,j)
                IsPCCell_candidate{j,1}(g)=k;
                g=g+1;
            else
            end
        else
        end
    end
end
%% check place fields and identify real PCs
PlaceFields=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
IsPCCell=cell(ExperimentInformation.Session,1);
for j=1:1:ExperimentInformation.Session
    g=1;
    for k=1:1:length(IsPCCell_candidate{j,1})
        PlaceFields{j,IsPCCell_candidate{j,1}(k)}=struct;
        [PlaceFields{j,IsPCCell_candidate{j,1}(k)}.fieldsMap,PlaceFields{j,IsPCCell_candidate{j,1}(k)}.fields]=analyses.placefield(ActivityMap{j,IsPCCell_candidate{j,1}(k)},'minBins',MinBins_field,'threshold',Cutoff_field,'minPeak',MinPeek_field,'binWidth',MapBinsize);
        if size(PlaceFields{j,IsPCCell_candidate{j,1}(k)}.fields,2)>=1
            for h=1:1:size(PlaceFields{j,IsPCCell_candidate{j,1}(k)}.fields,2)
                FieldSize(h)=PlaceFields{j,IsPCCell_candidate{j,1}(k)}.fields(h).area;
            end
            if max(FieldSize)<MaxFieldsize
                IsPCCell{j,1}(g)=IsPCCell_candidate{j,1}(k);
                g=g+1;
            else
            end
        else
        end
    end
end

%% Real PCs v.s. change
j=1; %session ID
TotalRun=1000;
TotalNumber=length(CellForAnalysis{1,1});
SelectCell=IsPCCell{session,1};
ChanceWin=zeros(TotalRun,1);
for w=1:1:TotalRun
    TotalNumber_PC_chance=1;
    while (TotalNumber_PC_chance<TotalNumber)
        for k=1:1:size(SelectCell,2)
            if TotalNumber_PC_chance<=TotalNumber
            i=SelectCell(k);
            SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
            SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
            SelectedFrame_firstHalf=SelectedFrame_filtered(find(SelectedFrame_filtered<=ExperimentInformation.FrameinEachSession/2));
            SelectedFrame_secondHalf=SelectedFrame_filtered(find(SelectedFrame_filtered>ExperimentInformation.FrameinEachSession/2));
            PositionTrain=NAT{1,j}(SelectedFrame_filtered,1:3);
            PositionTrain_firstHalf=NAT{1,j}(SelectedFrame_firstHalf,1:3);
            PositionTrain_secondHalf=NAT{1,j}(SelectedFrame_secondHalf,1:3);
            EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
            xmin=Shuffling_mininterval*ExperimentInformation.Trackingframerate;
            xmax=size(length(SelectedFrame_filtered),1)-xmin;
            ShiftFrame=round(xmin+rand(1,1)*(xmax-xmin));
            EventTrain_shuffled=circshift(EventTrain,ShiftFrame);
            EventTrain_shuffled_firstHalf=EventTrain_shuffled(1:size(PositionTrain_firstHalf,1),:);
            EventTrain_shuffled_secondHalf=EventTrain_shuffled(size(PositionTrain_firstHalf,1)+1:end,:);
            ActivityMap_shuffled=analyses.map(PositionTrain,EventTrain_shuffled,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',[-32 32 -32 32]);
            ActivityMapStats_shuffled=analyses.mapStatsPDF(ActivityMap_shuffled);
            SI_chance=ActivityMapStats_shuffled.content;
            Map_shuffled_firstHalf=analyses.map(PositionTrain_firstHalf,EventTrain_shuffled_firstHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime/2,'limits',[-32 32 -32 32]);
            Map_shuffled_secondHalf=analyses.map(PositionTrain_secondHalf,EventTrain_shuffled_secondHalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime/2,'limits',[-32 32 -32 32]);
            MapCorrelation_chance=analyses.spatialCrossCorrelation(Map_shuffled_firstHalf.z, Map_shuffled_secondHalf.z);            
            PlaceFields_chance=struct;
            [PlaceFields_chance.fieldsMap,PlaceFields_chance.fields]=analyses.placefield(ActivityMap_shuffled,'minBins',MinBins_field,'threshold',Cutoff_field,'minPeak',MinPeek_field,'binWidth',MapBinsize);
            FieldSize=[];
            if size(PlaceFields_chance.fields,2)>=1
                for h=1:1:size(PlaceFields_chance.fields,2)
                    FieldSize(h)=PlaceFields_chance.fields(h).area;
                end
            else
                FieldSize=0;
            end           
            Maxfield=max(FieldSize);            
            if SI_chance>SI_shuffled(i,Shuffling+1,j) && MapCorrelation_chance>MapCorrelation_shuffled(i,Shuffling+1,j) && Maxfield<MaxFieldsize
                ChanceWin(w,1)=ChanceWin(w,1)+1;
            else
            end
            TotalNumber_PC_chance=TotalNumber_PC_chance+1;
            disp([num2str(w),' ', num2str(TotalNumber_PC_chance)]);
            else
            end
        end
    end
end
%%
close all
figure
x0=100;
y0=100;
width=250;
height=200;
set(gcf,'position',[x0,y0,width,height])
histogram(ChanceWin/TotalNumber*100,10,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
hold on
plot([100*size(IsPCCell{j,1},2)/TotalNumber 100*size(IsPCCell{j,1},2)/TotalNumber],[0 0.3],'LineWidth',3,'color',[1 0 0]);
box off
set(gca, 'TickDir', 'out')
xlim([-0.5 50])
xticks([0:10:50])
ylim([0 0.4])
