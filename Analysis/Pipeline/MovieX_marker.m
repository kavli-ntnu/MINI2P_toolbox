%%
File_folder_raw=uigetdir; % munually select the raw file folder;

if File_folder_raw~=0
    
    disp('Data folder was founded');
else
    disp('Data folder does not exist');
end
cd (File_folder_raw); % the working folder will be changed to the data folder
%% find all tracking video
AVI_Files=dir([File_folder_raw,'\*.avi']);
[AVI_num ~]=size(AVI_Files);

if AVI_num~=0
    disp([num2str(AVI_num),' tracking video files were founded']);
    for i=1:1:AVI_num
        disp(AVI_Files(i).name);
    end
else
    disp('0 tracking videos files were not founded');
end

VideoPoor=cell(AVI_num,1);
k=1;
VideoSeq=[];
for i=1:1:AVI_num
    VideoPoor{i,1}=VideoReader(AVI_Files(i).name);
    for j=1:1:VideoPoor{i,1}.NumFrames
    VideoSeq(k,1)=i;
    VideoSeq(k,2)=j;
    k=k+1;
    end
end    


load('ExperimentInformation.mat')
load('GridCellAnalysis.mat')
load('NAAK.mat')

%%
TraceTime=50; %second
TraceFrame=round(TraceTime.*ExperimentInformation.FrameRate/2);
StartFrame=round(TraceFrame/2);
CMAP=[1 0 0;...
     0 1 0;
     1 0 1;
     0 1 1;
     0 0 1];
%% find all image files (GCamp)
P1_Files=dir([File_folder_raw,'\SingleSlice\P1\*.tif']);
[P1_num ~]=size(P1_Files);
P2_Files=dir([File_folder_raw,'\SingleSlice\P2\*.tif']);

%% Make Movie S10: Grid cells closed to each tother have different phases.

close all
aviobj2 = VideoWriter('MOVIE_S10.mp4','MPEG-4');
aviobj2.Quality=20;
aviobj2.FrameRate=75;
open(aviobj2)
%%
CellID=[113,146,128,127];
CellPL=CellID>ExperimentInformation.CellinEachPlan(1);

for Frame=1:1:25000

    figure (1)
    x0=0;
    y0=0;
    width=1280;
    height=720;
    set(gcf,'position',[x0,y0,width,height])     

    % plot full image
    for i=1:1:length(CellID)
        Xc(i)=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
        Yc(i)=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
       
    end
    C1X=round(mean(Xc));
    C1Y=round(mean(Yc));
    P2=imread([P2_Files(Frame).folder,'\',P2_Files(Frame).name]);
    subplot(5,8,[1:2,9:10],'align')
    imshow (P2);
    hold on
    plot([C1X-30 C1X+30 C1X+30 C1X-30 C1X-30] ,[C1Y-30 C1Y-30 C1Y+30 C1Y+30 C1Y-30],'Color',[1 1 1 0.5],'LineWidth',1);
    hold on
    plot([243 243],[7 7+103],'Color',[1 1 1],'LineWidth',2);
    hold on
    text(225,85,'200 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])    
    xlim([1,256])
    ylim([1,256])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
%     set(gca, 'XDir','reverse')
%     set(gca, 'YDir','reverse')
    title ('Full FOV (FOV2)','Color',[1 1 1]); 
    subplot(5,8,[3:4,11:12],'align')
    imshow (P2);
    hold on
    for i=1:1:length(CellID)
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            P=convhull(Xpoint,Ypoint); 
            plot(Xpoint(P),Ypoint(P),'color',[CMAP(i,:),0.8],'LineWidth',1.5);
            text(Xcenter,Ycenter,num2str(i),'Color', CMAP(i,:));
    end
    plot([C1X C1X]+27,[C1Y-12 C1Y+12]-15,'Color',[1 1 1],'LineWidth',2);
    hold on
    text(C1X+23,C1Y-9,'25 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])    
    xlim([C1X-30,C1X+30])
    ylim([C1Y-30,C1Y+30])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
%     set(gca, 'XDir','reverse')
%     set(gca, 'YDir','reverse')
    title ('Zoom in','Color',[1 1 1]);
% plot the spikes of select cells on the animal tracking vidio    
    subplot(5,8,[17:20,25:28,33:36],'align')
    Tracking=read(VideoPoor{VideoSeq(Frame*2,1),1},VideoSeq(Frame*2,2));
    imshow(Tracking);
    hold on
    for i=1:1:length(CellID)
        SelectSpike_full=NAAK{1,1}(:,4*CellID(i)+12)>0;
        Spike_normalized=NAAK{1,1}(SelectSpike_full,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                            
        SelectSpike_crop=SelectSpike_full(1:Frame*2);
        Spike_normalized_crop=NAAK{1,1}(SelectSpike_crop,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                            
        scatter(512+1/ExperimentInformation.BoxScale*NAAK{1,1}(SelectSpike_crop,2) ,512-1/ExperimentInformation.BoxScale*NAAK{1,1}(SelectSpike_crop,3),200.*Spike_normalized_crop,CMAP(i,:),'filled','LineWidth',0.5,'MarkerFaceAlpha',0.4);        
        hold on
    end
    plot([50 950],[1124 1124],'Color',[1 1 1],'LineWidth',3);
    hold on
    text(450,1070,'80 cm','Color',[1 1 1],'FontSize',12)
    daspect([1 1 1]); 
    xlim([1 1024])
    ylim([1 1125])
%         box on
    axis off
    hold off
    box off
    title (['TimeStamp: ',num2str(Frame/15,'%.2f'),' s'],'Color',[1 1 1]);    
    

% plot the spike trace for single cell
for i=1:1:length(CellID)
    subplot(5,8,5+8*(i-1),'align')
    SelectSpike_full=NAAK{1,1}(:,4*CellID(i)+12)>0;
    Spike_normalized=NAAK{1,1}(SelectSpike_full,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                            
    SelectSpike_crop=SelectSpike_full(1:Frame*2);
    Spike_normalized_crop=NAAK{1,1}(SelectSpike_crop,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                              
    plot(NAAK{1,1}(1:Frame*2,2),NAAK{1,1}(1: Frame*2,3),'color', [0.3 0.3 0.3],'LineWidth',0.2);
            set(gca, 'Color', 'None');
            axis off;
            daspect([1 1 1]);  
            hold on      
    scatter(NAAK{1,1}(SelectSpike_crop,2) ,NAAK{1,1}(SelectSpike_crop,3),30.*Spike_normalized_crop,CMAP(i,:),'filled','LineWidth',0.5,'MarkerFaceAlpha',0.9);        
    daspect([1 1 1]); 
    xlim([-40 40])
    ylim([-40 40])
    xticks([-40 0 40])
    yticks([-40 0 40])
    ax=gca;
    ax.XColor([1 1 1]);
    ax.YColor([1 1 1]);
%     axis equal
    box on
    axis on
    hold off
    title(['Cell ',num2str(i),' (ID:' num2str(CellID(i)),' )'],'Color',CMAP(i,:))    
%     title ('Cal.events plot','Color',CMAP(i,:));
end

    
% show the zoom-in image of the cell
    for i=1:1:length(CellID)
        subplot(5,8,6+8*(i-1),'align')
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            P=convhull(Xpoint,Ypoint); 
            imshow(P2);
            hold on
            plot(Xpoint(P),Ypoint(P),'color',[CMAP(i,:),0.8],'LineWidth',2);         
            hold off
            axis off
            xlim([Xcenter-15,Xcenter+15])
            ylim([Ycenter-15,Ycenter+15])
%             title ('Cal.events plot','Color',CMAP(i,:));
            title(['Cell ',num2str(i),' (ID:' num2str(CellID(i)),' )'],'Color',CMAP(i,:))    
            camroll(-90)
    end

% plot the cell calcium signal

for i=1:1:length(CellID)
    subplot(5,8,[7:8]+8*(i-1),'align')
    TraceRaw=NAAK{1,1}(:,4*CellID(i)+9);  
    SpikeRaw=NAAK{1,1}(:,4*CellID(i)+12); 
    Trace_Filtered=TraceRaw(~isnan(TraceRaw));
    Trace_Filtered=smooth(Trace_Filtered,2);
    Spike_Filtered=SpikeRaw(~isnan(SpikeRaw))>1;
%     Spike_Filtered=SpikeRaw(~isnan(SpikeRaw))./max(SpikeRaw(~isnan(SpikeRaw)));
    TimeStamp=NAAK{1,1}(~isnan(TraceRaw),1);
    Trace_Filtered_nor=(Trace_Filtered-min(Trace_Filtered))./(max(Trace_Filtered)-min(Trace_Filtered));
    if 1-StartFrame+Frame>=1
    plot(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame),Trace_Filtered_nor([1-StartFrame:TraceFrame-StartFrame]+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    hold on
    ylim([0 1.15])
%     plot(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame),1.05-0.2*Spike_Filtered([1-StartFrame:TraceFrame-StartFrame]+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    scatter(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame),1.2-0.2*Spike_Filtered([1-StartFrame:TraceFrame-StartFrame]+Frame),10,CMAP(i,:),'|','LineWidth',2);
    xlim([min(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame)) max(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame))])
    hold on
    else
    plot(TimeStamp(1:TraceFrame-StartFrame+Frame),Trace_Filtered_nor(1:TraceFrame-StartFrame+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    hold on
    ylim([0 1.15])
%     plot(TimeStamp(1:TraceFrame-StartFrame+Frame),1.1-0.2*Spike_Filtered(1:TraceFrame-StartFrame+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    scatter(TimeStamp(1:TraceFrame-StartFrame+Frame),1.2-0.2*Spike_Filtered(1:TraceFrame-StartFrame+Frame),10,CMAP(i,:),'|','LineWidth',2);
    xlim([max(TimeStamp(1:TraceFrame-StartFrame+Frame))-TraceTime max(TimeStamp(1:TraceFrame-StartFrame+Frame))])
    hold on    
    end          
    plot([TimeStamp(1+Frame+2),TimeStamp(1+Frame+2)],[0 1],'color',[1 1 1],'LineWidth',2);
    box off
    axis off
    hold off
    title('dF/F and calcium event location','Color',CMAP(i,:))   
end
%     set(gcf,'color',[0 0 0]);
%     set(gca,'color',[1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
CurrentFrame=getframe(gcf);
writeVideo(aviobj2,CurrentFrame);
Frame
end
close(aviobj2);
 
