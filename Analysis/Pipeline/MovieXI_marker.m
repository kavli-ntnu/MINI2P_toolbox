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
color_scheme_npg = [...
    0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
         0    0.6275    0.5294; ...
    0.4353    0.5294    0.7333; ...
    0.9529    0.6078    0.4980; ...
    0.5686    0.8196    0.7608; ...
    0.5176    0.5686    0.7059; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];
CMAP=color_scheme_npg;
%% find all image files (GCamp)
P1_Files=dir([File_folder_raw,'\SingleSlice\P1\*.tif']);
[P1_num ~]=size(P1_Files);
P2_Files=dir([File_folder_raw,'\SingleSlice\P2\*.tif']);

%% Make Movie S11: Grid cells with projection to hippocumpus (tdTomato labeling).


close all

aviobj3 = VideoWriter('MOVIE_S9.mp4','MPEG-4');
aviobj3.Quality=20;
aviobj3.FrameRate=75;
open(aviobj3)
%%
CellID=[39,34,30,118,129,136];;
CellPL=CellID>ExperimentInformation.CellinEachPlan(1);
%%
for Frame=1:1:25000

figure (1)
    x0=0;
    y0=0;
    width=1280;
    height=720;
    set(gcf,'position',[x0,y0,width,height])

    % plot full image
    subplot(6,10,[1:2,11:12,21:22],'align')
    P1=imread([P1_Files(Frame).folder,'\',P1_Files(Frame).name]);
    imshow (P1);
    hold on
    for i=1:1:length(CellID)
        if CellPL(i)==0
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            P=convhull(Xpoint,Ypoint); 
            plot([Xcenter-15 Xcenter+15 Xcenter+15 Xcenter-15 Xcenter-15] ,[Ycenter-15 Ycenter-15 Ycenter+15 Ycenter+15 Ycenter-15],'color',[CMAP(i,:),0.7],'LineWidth',1);  
            text(Xcenter-30,Ycenter+20,num2str(i),'Color', CMAP(i,:),'fontweight', 'bold' );
        else
        end
    end
    hold on
    plot([243 243],[7 7+103],'Color',[1 1 1],'LineWidth',2);
    hold on
    text(227,85,'200 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
%     set(gca, 'XDir','reverse')
%     set(gca, 'YDir','reverse')
    title ('P1: -140 micron (GCaMP6)','Color',[1 1 1]);
    
    % plot full image (Channel 2)

    subplot(6,10,[3:4,13:14,23:24],'align')
    P1C2=imread([P1C2_Files(1).folder,'\',P1C2_Files(1).name]);
    imshow (P1C2);
    hold on
    for i=1:1:length(CellID)
        if CellPL(i)==0
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            P=convhull(Xpoint,Ypoint); 
            plot([Xcenter-15 Xcenter+15 Xcenter+15 Xcenter-15 Xcenter-15] ,[Ycenter-15 Ycenter-15 Ycenter+15 Ycenter+15 Ycenter-15],'color',[CMAP(i,:),0.7],'LineWidth',1);  
            text(Xcenter-30,Ycenter+20,num2str(i),'Color', CMAP(i,:),'fontweight', 'bold' );
        else
        end
    end
    hold on
    plot([243 243],[7 7+103],'Color',[1 1 1],'LineWidth',2);
    hold on
    text(227,85,'200 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
%     set(gca, 'XDir','reverse')
%     set(gca, 'YDir','reverse')
    title ('P1: -140 micron (tdTomato)','Color',[255 189 000]/255);
    
    
    
    
    
    
    
    subplot(6,10,[31:32,41:42,51:52],'align')
    P2=imread([P2_Files(Frame).folder,'\',P2_Files(Frame).name]);
    imshow (P2);
    hold on
    for i=1:1:length(CellID)
        if CellPL(i)==1
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            P=convhull(Xpoint,Ypoint);             
            plot([Xcenter-15 Xcenter+15 Xcenter+15 Xcenter-15 Xcenter-15] ,[Ycenter-15 Ycenter-15 Ycenter+15 Ycenter+15 Ycenter-15],'color',[CMAP(i,:),0.7],'LineWidth',1);  
            text(Xcenter-30,Ycenter+20,num2str(i),'Color', CMAP(i,:),'fontweight', 'bold' );
        else
        end
    end
    hold on
    plot([243 243],[7 7+103],'Color',[1 1 1],'LineWidth',2);
    hold on
    text(227,85,'200 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
    title ('P2: -100 micron (GCaMP6)','Color',[1 1 1]);
  
    
    subplot(6,10,[33:34,43:44,53:54],'align')
    P2C2=imread([P2C2_Files(1).folder,'\',P2C2_Files(1).name]);
    imshow (P2C2);
    hold on
    for i=1:1:length(CellID)
        if CellPL(i)==1
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            P=convhull(Xpoint,Ypoint);             
            plot([Xcenter-15 Xcenter+15 Xcenter+15 Xcenter-15 Xcenter-15] ,[Ycenter-15 Ycenter-15 Ycenter+15 Ycenter+15 Ycenter-15],'color',[CMAP(i,:),0.7],'LineWidth',1);  
            text(Xcenter-30,Ycenter+20,num2str(i),'Color', CMAP(i,:),'fontweight', 'bold' );
        else
        end
    end
    hold on
    plot([243 243],[7 7+103],'Color',[1 1 1],'LineWidth',2);
    hold on
    text(227,85,'200 µm','Color',[1 1 1],'FontSize',10)
    hold off
    axis on
    xticks([])
    yticks([])
    camroll(-90)
    set(gca, 'Color',[1 1 1])
    title ('P2: -100 micron (tdTomato)','Color',[255 189 000]/255);
    
    
 for i=1:1:length(CellID)
    subplot(6,10,5+10*(i-1),'align')
    SelectSpike_full=NAAK{1,1}(:,4*CellID(i)+12)>0;
    Spike_normalized=NAAK{1,1}(SelectSpike_full,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                            
    SelectSpike_crop=SelectSpike_full(1:Frame*2);
    Spike_normalized_crop=NAAK{1,1}(SelectSpike_crop,CellID(i)*4+12)./max(NAAK{1,1}(SelectSpike_full,CellID(i)*4+12));                              
    plot(NAAK{1,1}(1:Frame*2,2),NAAK{1,1}(1: Frame*2,3),'color', [0.4 0.4 0.4],'LineWidth',0.2);
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
    title(['Cell ',num2str(i),' ( ID:' num2str(CellID(i)),' )'],'Color',CMAP(i,:))     
    
    
    
    
    
%   plot single cell activity (GCaMP6)  
     
    for i=1:1:length(CellID)
        subplot(6,10,6+10*(i-1),'align')
        if CellPL(i)==0
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            P=convhull(Xpoint,Ypoint); 
            imshow(P1);
            hold on
            plot(Xpoint(P),Ypoint(P),'color',[CMAP(i,:),0.8],'LineWidth',2);         
            hold off
            axis off
            xlim([Xcenter-15,Xcenter+15])
            ylim([Ycenter-15,Ycenter+15])
        else
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
        end
        title(['Cell ',num2str(i),' ( ID:' num2str(CellID(i)),' )'],'Color',CMAP(i,:))     
          camroll(-90)
    end   

    
%   plot single cell channel 2 average (tdTomata)  
     
    for i=1:1:length(CellID)
        subplot(6,10,7+10*(i-1),'align')
        if CellPL(i)==0
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            P=convhull(Xpoint,Ypoint); 
            imshow(P1C2);
            hold on
            plot(Xpoint(P),Ypoint(P),'color',[CMAP(i,:),0.8],'LineWidth',2);         
            hold off
            axis off
            xlim([Xcenter-15,Xcenter+15])
            ylim([Ycenter-15,Ycenter+15])
        else
            Xcenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(2))',256);
            Ycenter=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.med(1))',256);
            Xpoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.xpix)',256);
            Ypoint=mod(double(ExperimentInformation.CellStat{1,CellID(i)}.ypix)',256);
            P=convhull(Xpoint,Ypoint); 
            imshow(P2C2);
            hold on
            plot(Xpoint(P),Ypoint(P),'color',[CMAP(i,:),0.8],'LineWidth',2);         
            hold off
            axis off
            xlim([Xcenter-15,Xcenter+15])
            ylim([Ycenter-15,Ycenter+15])
        end
        title(['TdTomato (',num2str(redcell(CellID(i),2),'%.2f'),')'],'Color',CMAP(i,:))     
          camroll(-90)
    end      
    
%     
%     
% 
% 
% plot the cell calcium signal

for i=1:1:length(CellID)
    subplot(6,10,[8:9]+10*(i-1),'align')
    TraceRaw=NAAK{1,1}(:,4*CellID(i)+9);  
    SpikeRaw=NAAK{1,1}(:,4*CellID(i)+12); 
    Trace_Filtered=TraceRaw(~isnan(TraceRaw));
    Trace_Filtered=smooth(Trace_Filtered,3);
    Spike_Filtered=SpikeRaw(~isnan(SpikeRaw))>0.5;
    TimeStamp=NAAK{1,1}(~isnan(TraceRaw),1);
    Trace_Filtered_nor=(Trace_Filtered-min(Trace_Filtered))./(max(Trace_Filtered)-min(Trace_Filtered));
    if 1-StartFrame+Frame>=1
    plot(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame),Trace_Filtered_nor([1-StartFrame:TraceFrame-StartFrame]+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    hold on
    ylim([0 1.15])
    scatter(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame),1.2-0.2*Spike_Filtered([1-StartFrame:TraceFrame-StartFrame]+Frame),10,CMAP(i,:),'|','LineWidth',2);
    xlim([min(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame)) max(TimeStamp([1-StartFrame:TraceFrame-StartFrame]+Frame))])
    hold on
    else
    plot(TimeStamp(1:TraceFrame-StartFrame+Frame),Trace_Filtered_nor(1:TraceFrame-StartFrame+Frame),'color',CMAP(i,:),'LineWidth',1.5);
    hold on
    ylim([0 1.15])
    scatter(TimeStamp(1:TraceFrame-StartFrame+Frame),1.2-0.2*Spike_Filtered(1:TraceFrame-StartFrame+Frame),10,CMAP(i,:),'|','LineWidth',2);
    xlim([max(TimeStamp(1:TraceFrame-StartFrame+Frame))-TraceTime max(TimeStamp(1:TraceFrame-StartFrame+Frame))])
    hold on    
    end                  
    plot([TimeStamp(1+Frame+2),TimeStamp(1+Frame+2)],[0 1],'color',[1 1 1],'LineWidth',2);
    box off
    axis off
    hold off
    title('dF/F and calcium event location ','Color',CMAP(i,:))   
end



 end
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
CurrentFrame=getframe(gcf);
writeVideo(aviobj3,CurrentFrame);
Frame
end
close(aviobj3)

 
