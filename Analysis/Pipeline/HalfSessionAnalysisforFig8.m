close all;
%%
GridCellAnalysis.IsGridCell{1,1}=StichingPoor{2,1}.GridCellAnalysis.GridCellAnalysis.IsGridCell{1,1};
%%
% paremater settings
MapSmooth=2.5;
MapBinsize=2.5; % cm
SpeedThreadhold=2.5; %cm/mm
MinEventCount=100;
MinSNR=3;
Shuffling=1000;
Shuffling_mininterval=30;
MinTime=0.1;
CorrectionScan=[1 0 -1 -2];
Limit=[-40 40 -40 40];
Radii=[5 5];
GridBestShift=zeros(ExperimentInformation.Session,ExperimentInformation.TotalCell);
Gridscore_raw=zeros(length(size(GridCellAnalysis.IsGridCell{1,1},2)),3);
GridsStat_full=cell(length(size(GridCellAnalysis.IsGridCell{1,1},2)),1);
ActivityMap_firsthalf=cell(length(size(GridCellAnalysis.IsGridCell{1,1},2)),1);
ActivityMap_secondthalf=cell(length(size(GridCellAnalysis.IsGridCell{1,1},2)),1);
AutocorrelationMap_firsthalf=cell(length(size(GridCellAnalysis.IsGridCell{1,1},2)),1);
AutocorrelationMap_secondthalf=cell(length(size(GridCellAnalysis.IsGridCell{1,1},2)),1);
RadiRange=[5 5];

%% Colors used in Science journals
% Taken from https://ggsci.net/index.html 
color_scheme_aaas = [...
    0.2314    0.2863    0.5725; ...
    0.9333         0         0; ...
         0    0.5451    0.2706; ...
    0.3882    0.0941    0.4745; ...
         0    0.5098    0.5020; ...
    0.7333         0    0.1294; ...
    0.3725    0.3333    0.6078; ...
    0.6353         0    0.3373; ...
    0.5020    0.5059    0.5020; ...
    0.1059    0.0980    0.0980];


%%
j=1;
for k=1:1:size(GridCellAnalysis.IsGridCell{1,1},2)
    i=GridCellAnalysis.IsGridCell{1,1}(k);
    SelectedFrame_raw=find(~isnan(NAT{1,j}(:,4*i+10)));
    SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
    SelectedFrame_filtered_firsthalf=SelectedFrame_filtered(SelectedFrame_filtered<ExperimentInformation.FrameinEachSession/2);
    SelectedFrame_filtered_secondhalf=SelectedFrame_filtered(SelectedFrame_filtered>ExperimentInformation.FrameinEachSession/2);
    Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)>=0));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold   
    EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
    EventTrain_firsthalf=NAT{1,j}(SelectedFrame_filtered_firsthalf,[1 4*i+12]);
    EventTrain_secondhalf=NAT{1,j}(SelectedFrame_filtered_secondhalf,[1 4*i+12]);   
    PositionTrain_firsthalf=NAT{1,j}(SelectedFrame_filtered_firsthalf,1:3);
    PositionTrain_secondthalf=NAT{1,j}(SelectedFrame_filtered_secondhalf,1:3);
    DirectionTrain_firsthald=NAT{1,j}(SelectedFrame_filtered_firsthalf,4);
    DirectionTrain_secondthalf=NAT{1,j}(SelectedFrame_filtered_secondhalf,4);
    Shift=GridCellAnalysis.GridBestShift(i);
    
    PositionTrain_firsthalf(:,2)=PositionTrain_firsthalf(:,2)+Shift.* cos(DirectionTrain_firsthald * pi/180);
    PositionTrain_firsthalf(:,3)=PositionTrain_firsthalf(:,3)+Shift.* sin(DirectionTrain_firsthald * pi/180);
    
    PositionTrain_secondthalf(:,2)=PositionTrain_secondthalf(:,2)+Shift.* cos(DirectionTrain_secondthalf.* pi/180);
    PositionTrain_secondthalf(:,3)=PositionTrain_secondthalf(:,3)+Shift.* sin(DirectionTrain_secondthalf.* pi/180);    
    
    
    
    
    [Gridscore_full(k,1),GridsStat_full{k,1},GridCenter,~,GridsScoreRadius]=analyses.gridnessScore(GridCellAnalysis.AutocorrelationMap{1,i});  
    ActivityMap_firsthalf{k,1}=analyses.map(PositionTrain_firsthalf,EventTrain_firsthalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',0,'limits',Limit);                                 
    AutocorrelationMap_firsthalf{k,1}=analyses.autocorrelation(ActivityMap_firsthalf{k,1}.z);               
%     Gridscore_full(k,2)=BNTanalysis_modified.gridnessScoreinSession(AutocorrelationMap_firsthalf{k,1}, GridCenter, GridsScoreRadius, RadiRange);   
    Gridscore_full(k,2)=analyses.gridnessScore(AutocorrelationMap_firsthalf{k,1});  
    ActivityMap_secondthalf{k,1}=analyses.map(PositionTrain_secondthalf,EventTrain_secondhalf,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',0,'limits',Limit);                                 
    AutocorrelationMap_secondthalf{k,1}=analyses.autocorrelation(ActivityMap_secondthalf{k,1}.z);               
    Gridscore_full(k,3)=BNTanalysis_modified.gridnessScoreinSession(AutocorrelationMap_secondthalf{k,1}, GridCenter, GridsScoreRadius, RadiRange);
%     Gridscore_full(k,3)=analyses.gridnessScore(AutocorrelationMap_secondthalf{k,1});  

    k
end
    
%% 
[h,p_ttest] = ttest(Gridscore_full(:,2),Gridscore_full(:,3));   
%%
GridCell_ineachhalf=zeros(size(GridCellAnalysis.IsGridCell{1,1},2),2);
for i=1:1:size(GridCellAnalysis.IsGridCell{1,1},2)
    if Gridscore_full(i,2)>GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{1,1}(1,i),1001)
        GridCell_ineachhalf(i,1)=1;
    else
    end
    if Gridscore_full(i,3)>GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{1,1}(1,i),1001)
        GridCell_ineachhalf(i,2)=1;
    else
    end
end

Gridcell_still_firsthalf=find(GridCell_ineachhalf(:,1)==1);
Gridcell_still_secondhalf=find(GridCell_ineachhalf(:,2)==1);
Gridcell_onlyinfirsthalf=Gridcell_still_firsthalf(~ismember(Gridcell_still_firsthalf,Gridcell_still_secondhalf));
Gridcell_inbothhalf=Gridcell_still_firsthalf(ismember(Gridcell_still_firsthalf,Gridcell_still_secondhalf));
%%
GS_shuffled=GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{1,1},1:1000);
close all
figure
    x0=10;
    y0=10;
    width=250;
    height=250;
    set(gcf,'position',[x0,y0,width,height])
% subplot(1,3,1)
yyaxis left
histogram(GS_shuffled(:),'Normalization','probability','BinWidth',0.05,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
set(gca,  'XColor','k', 'YColor','k')
ylim([0 0.15])
yyaxis right
% histogram(Gridscore_full(:,1),'Normalization','probability','BinWidth',0.15,'DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 0 0] );
hold on
histogram(Gridscore_full(:,2),'Normalization','probability','BinWidth',0.2,'DisplayStyle','stairs','LineWidth',2,'EdgeColor',color_scheme_aaas(2,:));
hold on
histogram(Gridscore_full(:,3),'Normalization','probability','BinWidth',0.18,'DisplayStyle','stairs','LineWidth',2,'EdgeColor',color_scheme_aaas(1,:));

ylim([0 0.24])
xlim([-1.2,1.8])
xticks([-1:0.5:2])
box off 
set(gca, 'TickDir', 'out')

%%
figure
x0=100;
y0=100;
width=100;
height=250;
set(gcf,'position',[x0,y0,width,height])
bar(1,length(Gridcell_still_firsthalf)./size(GridCellAnalysis.IsGridCell{1,1},2),'FaceColor',color_scheme_aaas(2,:),'EdgeColor','none');
hold on
bar(2,length(Gridcell_still_secondhalf)./size(GridCellAnalysis.IsGridCell{1,1},2),'FaceColor',color_scheme_aaas(1,:),'EdgeColor','none');
box off 
set(gca, 'TickDir', 'out')
xticks([1 2]);
yticks([0:0.2:0.6]);
ylim([0 0.6]);
%%
selectcell=27;
CellID=Gridcell_inbothhalf(selectcell);
% CellID=Gridcell_still_firsthalf(selectcell);
j=1;
close all
figure

x0=100;
y0=100;
width=280;
height=600;
set(gcf,'position',[x0,y0,width,height])
CMP=WJplots.CMP.inferno(256);
    i=GridCellAnalysis.IsGridCell{1,1}(CellID);
    SelectedFrame_raw=find(~isnan(NAT{1,j}(:,4*i+10)));
    SelectedFrame_filtered=intersect(find(~isnan(NAT{1,j}(:,4*i+10))),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
    SelectedFrame_filtered_firsthalf=SelectedFrame_filtered(SelectedFrame_filtered<ExperimentInformation.FrameinEachSession/2);
    SelectedFrame_filtered_secondhalf=SelectedFrame_filtered(SelectedFrame_filtered>ExperimentInformation.FrameinEachSession/2);
    Event_filtered=intersect(find(NAT{1,j}(:,4*i+12)>0),find(NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold   
    EventTrain=NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
    EventTrain_firsthalf=NAT{1,j}(SelectedFrame_filtered_firsthalf,[1 4*i+12]);
    Event_filtered_firsthalf=EventTrain_firsthalf(EventTrain_firsthalf(:,2)>0,:);
    EventTrain_secondhalf=NAT{1,j}(SelectedFrame_filtered_secondhalf,[1 4*i+12]); 
    Event_filtered_secondhalf=EventTrain_secondhalf(EventTrain_secondhalf(:,2)>0,:);
    PositionTrain_firsthalf=NAT{1,j}(SelectedFrame_filtered_firsthalf,1:3);
    PositionTrain_secondhalf=NAT{1,j}(SelectedFrame_filtered_secondhalf,1:3);


subplot(4,2,1)

Max=max(Event_filtered_firsthalf(:,2));
plot(PositionTrain_firsthalf(:,2),PositionTrain_firsthalf(:,3),'color', [0.6 0.6 0.6],'LineWidth',0.5);
hold on
scatter(PositionTrain_firsthalf(EventTrain_firsthalf(:,2)>0,2),PositionTrain_firsthalf(EventTrain_firsthalf(:,2)>0,3),20*(Event_filtered_firsthalf(:,2)./Max),color_scheme_aaas(2,:),'filled','MarkerFaceAlpha',0.8)  
ylim([-41,41])
xlim([-41,41])
title([num2str(i),'- frist half']);
daspect([1 1 1]); 
box off
axis off

subplot(4,2,2)

Max=max(Event_filtered_secondhalf(:,2));
plot(PositionTrain_secondhalf(:,2),PositionTrain_secondhalf(:,3),'color', [0.6 0.6 0.6],'LineWidth',0.5);
hold on
scatter(PositionTrain_secondhalf(EventTrain_secondhalf(:,2)>0,2),PositionTrain_secondhalf(EventTrain_secondhalf(:,2)>0,3),20*(Event_filtered_secondhalf(:,2)./Max),color_scheme_aaas(2,:),'filled','MarkerFaceAlpha',0.8)  
ylim([-41,41])
xlim([-41,41])
title([num2str(i),'- second half']);
daspect([1 1 1]); 
box off
axis off



subplot(4,2,3)
MAP=ActivityMap_firsthalf{CellID,1}.z;
MAX=prctile(MAP(:),99.5);
imagesc(flipud(MAP./MAX),'AlphaData',flipud(MAP)>0);
colormap(CMP)
caxis([0 1] ); 
ylim([0 size(MAP,1)])
xlim([0 size(MAP,2)])
title(['P:',num2str(max(MAP(:)),'%.2f')]);
daspect([1 1 1]); 
box off
axis off

subplot(4,2,4)
MAP=ActivityMap_secondthalf{CellID,1}.z;
MAX=prctile(MAP(:),99.5);
imagesc(flipud(MAP./MAX),'AlphaData',flipud(MAP)>0);
colormap(CMP)
caxis([0 1] ); 
ylim([0 size(MAP,1)])
xlim([0 size(MAP,2)])
title(['P:',num2str(max(MAP(:)),'%.2f')]);
daspect([1 1 1]); 
box off
axis off

h=subplot(4,2,5);
MAP=AutocorrelationMap_firsthalf{CellID,1};
% MAX=prctile(MAP(:),99.5);
imagesc(flipud(MAP));
colormap(h,jet)
caxis([-1 1] ); 
ylim([0 size(MAP,1)])
xlim([0 size(MAP,2)])
title(['G:',num2str(Gridscore_full(CellID,2),'%.2f')]);
daspect([1 1 1]); 
box off
axis off

h=subplot(4,2,6);
MAP=AutocorrelationMap_secondthalf{CellID,1};
% MAX=prctile(MAP(:),99.5);
imagesc(flipud(MAP));
colormap(h,jet)
caxis([-1 1] ); 
ylim([0 size(MAP,1)])
xlim([0 size(MAP,2)])
title(['G:',num2str(Gridscore_full(CellID,3),'%.2f')]);
daspect([1 1 1]); 
box off
axis off    





h=subplot(4,2,7:8);
Shuffled_example=GridCellAnalysis.GridScore_shuffled(i,1:1000);
H=histogram(Shuffled_example,20,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
hold on
Cut=GridCellAnalysis.GridScore_shuffled(i,1001);
plot([Gridscore_full(CellID,2) Gridscore_full(CellID,2)],[0 1.3*max(H.Values)],'color',color_scheme_aaas(2,:),'LineWidth',3);
hold on
plot([Gridscore_full(CellID,3) Gridscore_full(CellID,3)],[0 1.3*max(H.Values)],'color',color_scheme_aaas(1,:),'LineWidth',3);
hold on
plot([Cut Cut],[0 1.3*max(H.Values)],'color',color_scheme_aaas(3,:),'LineWidth',3);
box off
ylim([0 ceil(max(H.Values)/0.05)*0.05])
set(gca, 'TickDir', 'out')
xlim([-1.2 1.2])
ylim([0 0.2])
% h=subplot(4,2,8);
% Shuffled_example=GridCellAnalysis.GridScore_shuffled(i,1:1000);
% H=histogram(Shuffled_example,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
% hold on
% Cut=GridCellAnalysis.GridScore_shuffled(i,1001);
% plot([Gridscore_full(CellID,3) Gridscore_full(CellID,3)],[0 max(H.Values)],'color',[1 0 0],'LineWidth',3);
% hold on
% plot([Cut Cut],[0 max(H.Values)],'color',[0 0.7 0.7],'LineWidth',3);
% box off
% ylim([0 ceil(max(H.Values)/0.05)*0.05])
% set(gca, 'TickDir', 'out')