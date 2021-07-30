close all;
% paremater settings
MapSmooth=2.5;
MapBinsize=2.5; % cm
SpeedThreadhold=2.5; %cm/mm
MinEventCount=100;
MinSNR=3;
Shuffling=1000;
Shuffling_mininterval=30; %second
MinTime=0.1; %second
CorrectionScan=[1 0 -1 -2];
Limit=[-41 41 -41 41];
Radii=[5 5];
GridBestShift=zeros(ExperimentInformation.Session,ExperimentInformation.TotalCell);
%% shuffling
SelectedFrame_raw=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell
SelectedFrame_filtered=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell, filtered by some cratirial
ActivityMap=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
AutocorrelationMap=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
GridsStat=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %%filtered the map with 1) speed threadhold, moving condition or other critirial.
CellForAnalysis=cell(ExperimentInformation.Session,1);
GridScore_shuffled=zeros(ExperimentInformation.TotalCell,Shuffling+3,ExperimentInformation.Session);
for j=1:1:ExperimentInformation.Session
    k=1;
    for i=1:1:ExperimentInformation.TotalCell
        % for i=1;
        SelectedFrame_raw=find(~isnan(NAT{1,j}(:,4*i+10)));
        SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadholds       
        if length(Event_filtered)>MinEventCount && ExperimentInformation.CellSNR(i)>3 && ~ismember(i,ExperimentInformation.RepeatCell)
            EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
            CellForAnalysis{j,1}(k)=i;
            Gridscore_try=zeros(length(CorrectionScan),1);
            GridsStat_try=cell(length(CorrectionScan),1);
            GridCenter_try=zeros(length(CorrectionScan),1);
            GridsScoreRadius_try=cell(length(CorrectionScan),1);
            for p=1:1:length(CorrectionScan)
                PositionTrain_raw=NAT{1,j}(:,1:3);
                PositionTrain_raw(:,2)=circshift(PositionTrain_raw(:,2),CorrectionScan(p));
                PositionTrain_raw(:,3)=circshift(PositionTrain_raw(:,3),CorrectionScan(p));
                PositionTrain_try=PositionTrain_raw(SelectedFrame_filtered,:)       ;
                ActivityMap_try=analyses.map(PositionTrain_try,EventTrain,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',Limit);
                AutocorrelationMap_try=analyses.autocorrelation(ActivityMap_try.z);
                
                [ Gridscore_try(p,1),GridsStat_try{p,1},GridCenter_try(p,1),~,GridsScoreRadius_try{p,1}]=analyses.gridnessScore(AutocorrelationMap_try);
            end
            if ~isempty(Gridscore_try)
                [Gridscore_best,BestShift]=max(Gridscore_try(:));
                %                 GridScore_shuffled(i,Shuffling+3,j)=Gridscore_best;
                GridsStat{j,i}=GridsStat_try{BestShift,1};
                Centerfield=GridCenter_try(BestShift,1);
                ScoreRadius=GridsScoreRadius_try{BestShift,1};
            end
            GridBestShift(j,i)=CorrectionScan(BestShift);
            PositionTrain_raw=NAT{1,j}(:,1:3);
            PositionTrain_raw(:,2)=circshift(PositionTrain_raw(:,2),CorrectionScan(BestShift));
            PositionTrain_raw(:,3)=circshift(PositionTrain_raw(:,3),CorrectionScan(BestShift));
            PositionTrain=PositionTrain_raw(SelectedFrame_filtered,:);
            ActivityMap{j,i}=analyses.map(PositionTrain,EventTrain,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',Limit);
            AutocorrelationMap{j,i}=analyses.autocorrelation(ActivityMap{j,i}.z);
            if ~isnan(Gridscore_best) &&  ~isempty(GridsStat{j,i}.spacing)
                Orientationcheck=GridsStat{j,i}.orientation;
                Distancecheck=GridsStat{j,i}.spacing;
                if abs(Orientationcheck(2)-Orientationcheck(1))>30 && abs(Orientationcheck(2)-Orientationcheck(1))<90 && abs(Orientationcheck(3)-Orientationcheck(2))>30 && abs(Orientationcheck(3)-Orientationcheck(2))<90
                    if max(Distancecheck)<2*min(Distancecheck)
                        GridScore_shuffled(i,Shuffling+3,j)=Gridscore_best;
                        for n=1:1:Shuffling
                            xmin=Shuffling_mininterval*ExperimentInformation.FrameRate;
                            xmax=size(length(SelectedFrame_filtered),1)-xmin;
                            ShiftFrame=round(xmin+rand(1,1)*(xmax-xmin));
                            EventTrain_shuffled=circshift(EventTrain,ShiftFrame);
                            EventTrain_shuffled(:,1)=EventTrain(:,1);
                            ActivityMap_shuffled=analyses.map(PositionTrain,EventTrain_shuffled,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',Limit);
                            AutocorrelationMap_shuffle=analyses.autocorrelation(ActivityMap_shuffled.z);
                            GridScore_shuffled(i,n,j)= analyses.gridnessScoreShuffled(AutocorrelationMap_shuffle, Centerfield, ScoreRadius, Radii);
                            %                                     GridScore_shuffled(i,n,j)= analyses.gridnessScore(AutocorrelationMap_shuffle);
                            disp([num2str(i),' ',num2str(n)]);
                        end
                        GridScore_shuffled(i,Shuffling+1,j)=prctile(GridScore_shuffled(i,1:Shuffling,j),95);
                        GridScore_shuffled(i,Shuffling+2,j)=prctile(GridScore_shuffled(i,1:Shuffling,j),99);
                    else
                    end
                else
                end
            else
            end
        else
        end
    end
end

%% identify Grid cells
IsGridCell=cell(ExperimentInformation.Session,1);
for j=1:1:ExperimentInformation.Session
    g=1;
    for i=1:1:ExperimentInformation.TotalCell
        if ~isempty(GridsStat{j,i})
            if GridScore_shuffled(i,Shuffling+3,j)>GridScore_shuffled(i,Shuffling+1,j)
                IsGridCell{j,1}(g)=i;
                g=g+1;
            else
            end
        else
        end
    end
end
%% save grid cell information
GridCellAnalysis=struct;
GridCellAnalysis.IsGridCell=IsGridCell;
GridCellAnalysis.GridBestShift=GridBestShift;
GridCellAnalysis.ActivityMap=ActivityMap;
GridCellAnalysis.AutocorrelationMap=AutocorrelationMap;
GridCellAnalysis.CellForAnalysis=CellForAnalysis;
GridCellAnalysis.GridScore_shuffled=GridScore_shuffled;
GridCellAnalysis.GridsStat=GridsStat;
GridCellAnalysis.Shuffling=Shuffling;
%% plot  whole-session calcium activity maps for all grid cells
close all
% CellPerPlot=50;
k=1;
Xrange=[0 500];
session=1;
figure
x0=10;
y0=10;
width=4000;
height=1400;
set(gcf,'position',[x0,y0,width,height])
j=1;
for i=1:1:size(GridCellAnalysis.IsGridCell{j,1},2)
    subplot(10,14,i,'align')
    MAP=GridCellAnalysis.ActivityMap{session,GridCellAnalysis.IsGridCell{j,1}(i)}.z;
    MAX=prctile(MAP(:),99);
    %             imagesc(flipud((MAP-min(min(MAP)))./(max(max(MAP))-min(min(MAP)))),'AlphaData',MAP>0);
    imagesc(flipud(MAP./MAX),'AlphaData',flipud(MAP./MAX));
    
    CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    colormap(parula)
    %             colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP,1)])
    xlim([0 size(MAP,2)])
    title(['#',num2str(GridCellAnalysis.IsGridCell{j,1}(i)),' GC:',num2str(GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{j,1}(i),Shuffling+3,j),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    daspect([1 1 1]);
    box off
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);

%% plot  whole-session autocorrelation maps of calcium activity for all grid celss
close all
% CellPerPlot=50;
k=1;
Xrange=[0 500];
session=1;
figure
x0=10;
y0=10;
width=4000;
height=1400;
set(gcf,'position',[x0,y0,width,height])
j=1;
for i=1:1:size(GridCellAnalysis.IsGridCell{j,1},2)
    subplot(10,14,i,'align')
    MAP=GridCellAnalysis.AutocorrelationMap{session,GridCellAnalysis.IsGridCell{j,1}(i)};
    imagesc(flipud((MAP+1)./2));
    CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP,1)])
    xlim([0 size(MAP,2)])
    title(['#',num2str(GridCellAnalysis.IsGridCell{j,1}(i)),' GC:',num2str(GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{j,1}(i),Shuffling+3,j),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    daspect([1 1 1]);
    box off
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);
%% plot  whole-session activity may and autocorrelation maps of calcium activity for all grid celss
close all
% CellPerPlot=50;
k=1;
Xrange=[0 500];
session=1;
figure
x0=10;
y0=10;
width=3440;
height=1200;
set(gcf,'position',[x0,y0,width,height])
j=1;
for i=1:1:size(GridCellAnalysis.IsGridCell{j,1},2)
    subplot(10,ceil(size(GridCellAnalysis.IsGridCell{j,1},2)/10),i,'align')
    MAP1=GridCellAnalysis.ActivityMap{session,GridCellAnalysis.IsGridCell{j,1}(i)}.z;
    MAP1(isnan(MAP1))=0;
    MAP1=MAP1./prctile(MAP1(:),100);
    MAP2=GridCellAnalysis.AutocorrelationMap{session,GridCellAnalysis.IsGridCell{j,1}(i)};
    MAP2=(MAP2-min(MAP2(:)))./(1-min(MAP2(:)));
    MAP1=imresize(MAP1,size(MAP2));
    MAP_combine=[MAP1,MAP2];
    
    imagesc(flipud(MAP_combine));
    CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP_combine,1)])
    xlim([0 size(MAP_combine,2)])
    title(['#',num2str(GridCellAnalysis.IsGridCell{j,1}(i)),' GC:',num2str(GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{j,1}(i),Shuffling+3,j),'%.2f'),' P:',num2str(max(max(MAP1)),'%.2f')]);
    daspect([1 1 1]);
    box off
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);

%%
Gridcell_overlapped=[159 217 280 301 317];
GridCellAnalysis.Gridcell_overlapped=Gridcell_overlapped;
%%
save ([ExperimentInformation.RawDataAddress,'\GridCellAnalysis.mat'],'GridCellAnalysis','-v7.3');
disp('Grid cell analysis was done!');