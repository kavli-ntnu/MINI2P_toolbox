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
Shiftcorrection=0;
SelectedFrame_raw=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell
SelectedFrame_filtered=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell); %% which frame in each session should be used to do analyses for each cell, filtered by some cratirial
ActivityMap=cell(ExperimentInformation.Session,ExperimentInformation.TotalCell);
for j=1:1:ExperimentInformation.Session
    for i=1:1:ExperimentInformation.TotalCell     
        SelectedFrame_raw=find(~isnan(NAAK{1,j}(:,4*i+10)));
        SelectedFrame_filtered=intersect(find(~isnan(NAAK{1,j}(:,4*i+10))),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
        SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAAK{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        Event_filtered=intersect(find(NAAK{1,j}(:,4*i+12)>0),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(NAAK{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold                                         
        PositionTrain=NAAK{1,j}(SelectedFrame_filtered,1:3);   
        PositionTrain(:,2)=circshift(PositionTrain(:,2),Shiftcorrection);
        PositionTrain(:,3)=circshift(PositionTrain(:,3),Shiftcorrection);
        EventTrain=NAAK{1,j}(SelectedFrame_filtered,[1 4*i+12]);
        if ~isempty(EventTrain)
        ActivityMap{j,i}=analyses.map(PositionTrain,EventTrain,'smooth',MapSmooth,'binWidth',MapBinsize,'minTime',MinTime,'limits',[-32 32 -32 32]);
        else
        end
    end
end
% plot  whole-session maps of calcium activity for all cells
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
for i=1:1:ExperimentInformation.TotalCell     
        subplot(10,29,i,'align')
        if ~isempty(ActivityMap{session,i})
            MAP=ActivityMap{session,i}.z;
            imagesc(flipud((MAP-min(min(MAP)))./(max(max(MAP))-min(min(MAP)))),'AlphaData',MAP>0);
            CMP=WJplots.CMP.inferno(256);
            colormap(CMP)
            caxis([0 1] ); 
            ylim([0 size(MAP,1)])
            xlim([0 size(MAP,2)])
            title(['#',num2str(i),' P:',num2str(max(max(MAP)),'%.2f')]);
            daspect([1 1 1]); 
            box off
            axis off
        else
        end

end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);