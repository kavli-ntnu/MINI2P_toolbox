close all
clear all
%% 
color_scheme_npg = [...
    0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
         0    0.6275    0.5294; ...
    0.2353    0.3294    0.5333; ...
    0.9529    0.6078    0.4980; ...
    0.5176    0.5686    0.7059; ...
    0.5686    0.8196    0.7608; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];

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

%% step 1 load the averaged image

load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.mat');
FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.csv');                     
Image_mean=ExperimentInformation.ImagingOptions.meanImg(:,1:256);
Image_mean_corrected=Image_mean;
% Image_mean_corrected=imwarp(Image_mean,TransformMatrix,'OutputView',imref2d(size(Image_mean)));  
h=imagesc((Image_mean_corrected-min(min(Image_mean)))/(1*prctile(Image_mean(:),99)-min(min(Image_mean))),[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
set(h,'AlphaData', Image_mean_corrected>0);
hold on
colormap(gca,gray)
daspect([1 1 1]);
xlabel('x')
ylabel('y')
camroll(-90)
set(gca, 'TickDir', 'out')   

%% plot all cells with PC colored
close all
CMP=jet(256);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
% for k=1:1:size(GridCellAnalysis.IsGridCell{1,1},2)
%     i=GridCellAnalysis.IsGridCell{1,1}(1,k);
for k=1:1:size(IsPCCell{1,1},2)
    i=IsPCCell{1,1}(1,k);
%     if i>ExperimentInformation.CellinEachPlan(1)
        Xpoint=mod(double(ExperimentInformation.CellStat{1,i}.xpix)',256);
        Ypoint=mod(double(ExperimentInformation.CellStat{1,i}.ypix)',256);
        Xpoint(Xpoint==0)=256;
        Ypoint(Ypoint==0)=256;

        ROI_pannel=zeros(size(Image_mean));
        for m=1:1:length(Xpoint)
            ROI_pannel(Ypoint(m),Xpoint(m))=1;
        end
%         ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));        
        ROI_pannel_corrected=ROI_pannel;             

        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);              
        P=convhull(Xpoint,Ypoint);
        CenterX=mean(Xpoint);
        CenterY=mean(Ypoint);
        fill(Xpoint(P),Ypoint(P),[0.8 0 0],'LineStyle','none','facealpha',.8); 
        text(CenterX,CenterY,num2str(i));
        hold on;    
%     else
%     end
end
   hold off
    box off
%     axis off
    axis square
    daspect([1 1 1]);
    xlim([0 size(Image_mean,1)]);
    xlabel('x')
    ylim([0 size(Image_mean,2)]);
    ylabel('y')
%     colormap(a3,gray)
    caxis([0 1] );
    camroll(180)
    set(gca,'YDir','reverse');
%     set(gcf,'color',[0 0 0]);
%     set(gca,'color',[0 0 0]);
    set(a3, 'XAxisLocation', 'top')
    set(a3, 'TickDir', 'out')
    xticks([])
    yticks([])
%% find bad frames

Xpostion=NAAK{1,1}(:,2);
Ypostion=NAAK{1,1}(:,3);
Distance= (Xpostion(2:end)-Xpostion(1:end-1)).^2+(Ypostion(2:end)-Ypostion(1:end-1)).^2;
Badframe=find(Distance>10);

%%
close all

% SeceletCell=[113,125,127,146,128];
% SeceletCell=[125,113,146,128,127];
% SeceletCell=[113,146,128,127];
SeceletCell=[438 290];
% SeceletCell=[187 40 72];
% SeceletCell=[307 373 311];
NumCell=length(SeceletCell);
CMP=jet(NumCell);
% CMP=color_scheme_aaas;
figure
width=200*NumCell;
height=500;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
j=1;
CMP=[1 0 1;...
     0 1 1;
     0 1 0;
     0 0 1;
     0 0 1]
     
     
    
for k=1:1:NumCell
    i=SeceletCell(k);    
    subplot(3,NumCell+1,k)
    AnimaPosition=NAAK{1,j}(:,2:3);
    AnimaPosition(Badframe,1)=nan;
    AnimaPosition(Badframe,2)=nan;
    SelectedFrame_filtered=intersect(find(~isnan(NAAK{1,j}(:,4*i+10))),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(NAAK{1,j}(:,4*i+12)>0),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    AnimaPosition=AnimaPosition(SelectedFrame_filtered,:);
    
    
    
    
    EventPosition=NAAK{1,j}(Event_filtered,2:3);
    Event=NAAK{1,j}(Event_filtered,4*i+12);
    Max=max(Event(:));
    plot(smooth(AnimaPosition(:,1),'lowess',5),smooth(AnimaPosition(:,2),'lowess',5),'color', [0.5 0.5 0.5],'LineWidth',0.3);
    hold on
    scatter(EventPosition(:,1),EventPosition(:,2),25*(Event./Max),CMP(k,:),'filled','MarkerFaceAlpha',0.6)  
    ylim([-41,41])
    xlim([-41,41])
    title(['Grid cell ',num2str(k)]);
    daspect([1 1 1]); 
    box off
    axis off
end   
       
subplot(3,NumCell+1,NumCell+1) 
plot(AnimaPosition(:,1),AnimaPosition(:,2),'color', [0.3 0.3 0.3],'LineWidth',0.5);
for k=NumCell:-1:1
    i=SeceletCell(k);
    AnimaPosition=NAAK{1,j}(:,2:3);
    AnimaPosition(Badframe,1)=nan;
    AnimaPosition(Badframe,2)=nan;
    SelectedFrame_filtered=intersect(find(~isnan(NAAK{1,j}(:,4*i+10))),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(NAAK{1,j}(:,4*i+12)>0),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    AnimaPosition=AnimaPosition(SelectedFrame_filtered,:);
    EventPosition=NAAK{1,j}(Event_filtered,2:3);
    Event=NAAK{1,j}(Event_filtered,4*i+12);
    Max=max(Event(:));
    hold on
    scatter(EventPosition(:,1),EventPosition(:,2),25*(Event./Max),CMP(k,:),'filled','MarkerFaceAlpha',0.6)  
    ylim([-41,41])
    xlim([-41,41])
    title(['Overlap']);
    daspect([1 1 1]); 
    box off
    axis off
end

MAP={};
MAP_normalized={};
for k=1:1:NumCell
    i=SeceletCell(k);    
    subplot(3,NumCell+1,k+NumCell+1)
    MAP{1,k}=ActivityMap{1,i}.z;
%     MAX=prctile(MAP{1,k}(:),100);
    MAX=max(MAP{1,k}(:));
%     MIN=prctile(MAP{1,k}(:),0);
    MIN=min(MAP{1,k}(:));
    MAP_normalized{1,k}=(MAP{1,k}-MIN)./(MAX-MIN);
%     MAP_normalized{1,k}(MAP_normalized{1,k}>1)=1;
%     MAP_normalized{1,k}(MAP_normalized{1,k}<0)=0;       
    h=imagesc(flipud(MAP_normalized{1,k}),'AlphaData',MAP_normalized{1,k}>0);   
    colormap(gca,WJplots.SingleColorCMP(CMP(k,:),256));
    rgbImage{1,k} = ind2rgb(round((MAP{1,k}-MIN)./(MAX-MIN).*255), WJplots.SingleColorCMP(CMP(k,:),256));
    caxis([0 1] ); 
    ylim([2 size(MAP{1,k},1)-2])
    xlim([2 size(MAP{1,k},2)-2])
% title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
    daspect([1 1 1]); 
    box off
    axis off
%     imwrite(MAP_normalized{1,k},'map3.tiff','WriteMode','append');
end

subplot(3,NumCell+1,2*(NumCell+1))        

Map_merge=zeros(size(rgbImage{1,1}));
for i=1:1:NumCell
    Map_merge=Map_merge+rgbImage{1,i};
end
imagesc(flipud(Map_merge));

    ylim([2 size(Map_merge,1)-2])
    xlim([2 size(Map_merge,2)-2])
    daspect([1 1 1]); 
    box off
    axis off
    
    
subplot(3,NumCell+1,2*(NumCell+1)+1:3*(NumCell+1))     
m=1;
Offset=-1;
    for k=1:1:NumCell
        i=k;
        plot(NeuronActiveMatrix.F_raw_IsCell(:,1),smooth((m-1)*Offset+NeuronActiveMatrix.detaF_F(:,i+2),3),'color',CMP(k,:),'Linewidth',1);
        hold on
        m=m+1;
    end
ylim([NumCell*Offset 2])
xlim([0,400])
box off
axis off
%% draw the outline of selected cells on averaged image

close all
figure
width=500;
height=500;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
ImageProjection=imread('Image_projection-4.tif');
% ImageProjection_corrected=imwarp(ImageProjection,TransformMatrix,'OutputView',imref2d(size(ImageProjection)));
ImageProjection_corrected=ImageProjection;
h=imagesc(double(ImageProjection_corrected)/256,[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
set(h,'AlphaData', 1);
colormap(gca,gray)
%     caxis([0 1] );
hold on
for k=1:1:NumCell
        i=SeceletCell(k);
        Xpoint=mod(double(ExperimentInformation.CellStat{1,i}.xpix)',256);
        Ypoint=mod(double(ExperimentInformation.CellStat{1,i}.ypix)',256);
        Xpoint(Xpoint==0)=256;
        Ypoint(Ypoint==0)=256;

        ROI_pannel=zeros(size(Image_mean));
        for m=1:1:length(Xpoint)
            ROI_pannel(Ypoint(m),Xpoint(m))=1;
        end
%         ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));  
        ROI_pannel_corrected=ROI_pannel;        

        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);              
        P=convhull(Xpoint,Ypoint);
        CenterX=mean(Xpoint);
        CenterY=mean(Ypoint);
        plot(Xpoint(P),Ypoint(P),'color',[CMP(k,:),0.8],'LineWidth',1.5);
        hold on
end
   hold off
    box off
%     axis off
    axis square
    daspect([1 1 1]);
    xlim([0 size(ImageProjection_corrected,1)]);
    xlabel('x')
    ylim([0 size(ImageProjection_corrected,2)]);
    ylabel('y')
% 
%     camroll(-90)
%     set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
    xticks([])
    yticks([])

