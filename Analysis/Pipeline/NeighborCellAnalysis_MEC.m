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

%% plot all cells with grid cells color 
close all
CMP=jet(256);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
for k=1:1:size(GridCellAnalysis.IsGridCell{1,1},2)
    i=GridCellAnalysis.IsGridCell{1,1}(1,k);
    if i>ExperimentInformation.CellinEachPlan(1)
        Xpoint=mod(double(ExperimentInformation.CellStat{1,i}.xpix)',256);
        Ypoint=mod(double(ExperimentInformation.CellStat{1,i}.ypix)',256);
        Xpoint(Xpoint==0)=256;
        Ypoint(Ypoint==0)=256;

        ROI_pannel=zeros(size(Image_mean));
        for m=1:1:length(Xpoint)
            ROI_pannel(Ypoint(m),Xpoint(m))=1;
        end
        ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));        
%         ROI_pannel_corrected=ROI_pannel;             

        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);              
        P=convhull(Xpoint,Ypoint);
        CenterX=mean(Xpoint);
        CenterY=mean(Ypoint);
        fill(Xpoint(P),Ypoint(P),[0.8 0 0],'LineStyle','none','facealpha',.8); 
        text(CenterX,CenterY,num2str(i));
        hold on;    
    else
    end
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
    camroll(-90)
    set(gca,'YDir','reverse');
%     set(gcf,'color',[0 0 0]);
%     set(gca,'color',[0 0 0]);
    set(a3, 'XAxisLocation', 'top')
    set(a3, 'TickDir', 'out')
    xticks([])
    yticks([])

%%
close all

% SeceletCell=[113,125,127,146,128];
% SeceletCell=[125,113,146,128,127];
SeceletCell=[113,128,127];
% SeceletCell=[367,288,324,54];
NumCell=length(SeceletCell);
CMP=jet(NumCell);
% CMP=color_scheme_aaas;
figure
width=200*NumCell;
height=600;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
j=1;
CMP=[1 1 0;...
     1 0 1;
     0 1 1;
     1 0 1;
     0 0 1]
     
     
    
for k=1:1:NumCell
    i=SeceletCell(k);    
    subplot(4,NumCell+1,k)
    Bestshift=GridCellAnalysis.GridBestShift(i);
    SelectedFrame_filtered=intersect(find(~isnan(NAAK{1,j}(:,4*i+10))),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(NAAK{1,j}(:,4*i+12)>0),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    AnimaPosition=NAAK{1,j}(SelectedFrame_filtered,2:3);
    EventPosition=NAAK{1,j}(Event_filtered,2:3);
    EventPosition(:,1)=circshift(EventPosition(:,1),Bestshift);
    EventPosition(:,2)=circshift(EventPosition(:,2),Bestshift);
    Event=NAAK{1,j}(Event_filtered,4*i+12);
    Max=max(Event(:));
    plot(AnimaPosition(:,1),AnimaPosition(:,2),'color', [0.3 0.3 0.3],'LineWidth',0.3);
    hold on
    scatter(EventPosition(:,1),EventPosition(:,2),15*(Event./Max),CMP(k,:),'filled','MarkerFaceAlpha',0.9)  
    ylim([-41,41])
    xlim([-41,41])
    title(['Grid cell ',num2str(k)]);
    daspect([1 1 1]); 
    box off
    axis off
end   
       
subplot(4,NumCell+1,NumCell+1) 
plot(AnimaPosition(:,1),AnimaPosition(:,2),'color', [0.3 0.3 0.3],'LineWidth',0.5);
for k=1:1:NumCell
    i=SeceletCell(k);
    Bestshift=GridCellAnalysis.GridBestShift(i);
    SelectedFrame_filtered=intersect(find(~isnan(NAAK{1,j}(:,4*i+10))),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(NAAK{1,j}(:,4*i+12)>0),find(NAAK{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(NAAK{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    AnimaPosition=NAAK{1,j}(SelectedFrame_filtered,2:3);
    EventPosition=NAAK{1,j}(Event_filtered,2:3);
    EventPosition(:,1)=circshift(EventPosition(:,1),Bestshift);
    EventPosition(:,2)=circshift(EventPosition(:,2),Bestshift);
    Event=NAAK{1,j}(Event_filtered,4*i+12);
    Max=max(Event(:));

    hold on
    scatter(EventPosition(:,1),EventPosition(:,2),15*(Event./Max),CMP(k,:),'filled','MarkerFaceAlpha',0.9)  
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
    subplot(4,NumCell+1,k+NumCell+1)
    MAP{1,k}=GridCellAnalysis.ActivityMap{1,i}.z;
%     MAX=prctile(MAP{1,k}(:),100);
    MAX=max(MAP{1,k}(:));
%     MIN=prctile(MAP{1,k}(:),0);
    MIN=min(MAP{1,k}(:));
    MAP_normalized{1,k}=(MAP{1,k}-MIN)./(MAX-MIN);
%     MAP_normalized{1,k}(MAP_normalized{1,k}>1)=1;
%     MAP_normalized{1,k}(MAP_normalized{1,k}<0)=0;       
    h=imagesc(flipud(MAP_normalized{1,k}),'AlphaData',flipud(MAP_normalized{1,k})>0);   
    colormap(gca,WJplots.SingleColorCMP(CMP(k,:),256));
    rgbImage{1,k} = ind2rgb(round((MAP{1,k}-MIN)./(MAX-MIN).*255), WJplots.SingleColorCMP(CMP(k,:),256));
    caxis([0 1] ); 
    ylim([0 size(MAP{1,k},1)])
    xlim([0 size(MAP{1,k},2)])
% title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
    daspect([1 1 1]); 
    box off
    axis off
    imwrite(MAP_normalized{1,k},'map.tiff','WriteMode','append');
end

subplot(4,NumCell+1,2*(NumCell+1))        

Map_merge=zeros(size(rgbImage{1,1}));
for i=1:1:NumCell
    Map_merge=Map_merge+rgbImage{1,i};
end
imagesc(flipud(Map_merge));

    ylim([0 size(Map_merge,1)])
    xlim([0 size(Map_merge,2)])
    daspect([1 1 1]); 
    box off
    axis off
    
    

    
AutoMAP={};
AutoMAP_normalized={};
for k=1:1:NumCell
    i=SeceletCell(k);    
    subplot(4,NumCell+1,k+2*(NumCell+1))
    AutoMAP{1,k}=GridCellAnalysis.AutocorrelationMap{1,i};
    MIN=min(AutoMAP{1,k}(:));
    AutoMAP_normalized{1,k}=(AutoMAP{1,k}-MIN)./(1-MIN);
%     MAP_normalized{1,k}(MAP_normalized{1,k}>1)=1;
%     MAP_normalized{1,k}(MAP_normalized{1,k}<0)=0;       
    h=imagesc(flipud(AutoMAP_normalized{1,k}),'AlphaData',flipud(AutoMAP_normalized{1,k})>0);   
    colormap(gca,WJplots.SingleColorCMP(CMP(k,:),256));
    rgbImage{1,k} = ind2rgb(round((AutoMAP{1,k}-MIN)./(1-MIN).*255), WJplots.SingleColorCMP(CMP(k,:),256));
    caxis([0 1] ); 
    ylim([0 size(AutoMAP{1,k},1)])
    xlim([0 size(AutoMAP{1,k},2)])
% title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
    daspect([1 1 1]); 
    box off
    axis off
    imwrite(AutoMAP_normalized{1,k},'Automap.tiff','WriteMode','append');
end

subplot(4,NumCell+1,3*(NumCell+1))        

Map_merge=zeros(size(rgbImage{1,1}));
for i=1:1:NumCell
    Map_merge=Map_merge+rgbImage{1,i};
end
imagesc(flipud(Map_merge));

    ylim([0 size(Map_merge,1)])
    xlim([0 size(Map_merge,2)])
    daspect([1 1 1]); 
    box off
    axis off    
    
    
    
    
    
    
    
subplot(4,NumCell+1,3*(NumCell+1)+1:4*(NumCell+1))     
m=1;
Offset=-4;
    for k=1:1:NumCell
        i=k;
        plot(NeuronActiveMatrix.F_raw_IsCell(:,1),smooth((m-1)*Offset+NeuronActiveMatrix.detaF_F(:,i+2),3),'color',CMP(k,:),'Linewidth',1);
        hold on
        m=m+1;
    end
ylim([NumCell*Offset 5])
xlim([0,100])
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
ImageProjection=imread('Image_projection.tif');
ImageProjection_corrected=imwarp(ImageProjection,TransformMatrix,'OutputView',imref2d(size(ImageProjection)));
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
        ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));             
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

    camroll(-90)
    set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
    xticks([])
    yticks([])
%% calculate the crosscorrelation

MapCrossCorrelation=cell(NumCell,NumCell);
for i=1:NumCell
    for j=1:NumCell
        MapCrossCorrelation{i,j}=preprocessing.crosscorrelation(GridCellAnalysis.ActivityMap{1,i}.z,GridCellAnalysis.ActivityMap{1,j}.z);
        if ~(i==j)
        imwrite(MapCrossCorrelation{i,j},'MapCrossCorrelation.tiff','WriteMode','append');
        else
        end
    end
end

%%
MapCrossCorrelation_averaged=zeros(size(MapCrossCorrelation{1,1}));

m=1;
for i=1:NumCell
    for j=1:NumCell
        if ~(i==j)         
        MapCrossCorrelation_averaged=MapCrossCorrelation_averaged+MapCrossCorrelation{i,j};
        m=m+1
        else
        end
    end
end
MapCrossCorrelation_averaged=MapCrossCorrelation_averaged./(m-1);
%%

CMP=WJplots.fire(256);
figure
width=500;
height=500;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
close all

h=imagesc(MapCrossCorrelation_averaged,[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
MIN=min(MapCrossCorrelation_averaged(:))
MAX=max(MapCrossCorrelation_averaged(:))
set(h,'AlphaData', 1);
colormap(gca,CMP)
caxis([-0.12 0.12] );colorbar
daspect([1 1 1]);
xlim([30-20,30+20]);
ylim([30-20,30+20]);
    set(gca, 'TickDir', 'out')
%%
OutCircleRadi=20;
Center=ceil(size(MapCrossCorrelation_averaged)/2);
CircleVector=[0:1:OutCircleRadi];
CrossCorrelation_Peak=zeros(m-1,OutCircleRadi+1);
AutoCorrelation_Peak=zeros(5,OutCircleRadi+1);
Distnace=zeros(size(MapCrossCorrelation{1,1}));
for i=1:1:size(Distnace,1)
    for j=1:1:size(Distnace,2)
        Distnace(i,j)=round(sqrt((i-Center(1))^2+(j-Center(2))^2));
    end
end




%%
p=1;
n=1;
for i=1:1:NumCell
    for j=1:1:NumCell
        if ~(i==j)   
        MAP=MapCrossCorrelation{i,j};
        MAP=MAP(:);
        for k=CircleVector
            A=find(Distnace==k);
            CrossCorrelation_Peak(p,k+1)=max(MAP(A));
        end
            p=p+1;
        else
            MAP=MapCrossCorrelation{i,j};
            MAP=MAP(:);
            for k=CircleVector
            A=find(Distnace==k);
            AutoCorrelation_Peak(n,k+1)=max(MAP(A));
            end
            n=n+1;
        end
    end
end
        
%%
for i=1:1:size(AutoCorrelation_Peak,1)
    AutoCorrelation_Peak(i,:)=(AutoCorrelation_Peak(i,:)-min(AutoCorrelation_Peak(i,:)))./(max(AutoCorrelation_Peak(i,:))-min(AutoCorrelation_Peak(i,:)));
end

for i=1:1:size(CrossCorrelation_Peak,1)
    CrossCorrelation_Peak(i,:)=(CrossCorrelation_Peak(i,:)-min(CrossCorrelation_Peak(i,:)))./(max(CrossCorrelation_Peak(i,:))-min(CrossCorrelation_Peak(i,:)));
end  
    
    
AutoCorrelation_mean=mean(AutoCorrelation_Peak,1);
CrossCorrelation_mean=mean(CrossCorrelation_Peak,1);

AutoCorrelation_std=std(AutoCorrelation_Peak,1);
CrossCorrelation_std=std(CrossCorrelation_Peak,1);

close all
figure
width=275;
height=175;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])

plot(CircleVector,AutoCorrelation_mean,'color',[0 0 1],'LineWidth',1);
hold on

x2 = [CircleVector, fliplr(CircleVector)];
inBetween = [AutoCorrelation_mean-AutoCorrelation_std, fliplr(AutoCorrelation_mean+AutoCorrelation_std)];
fill(x2, inBetween, [0 0 1],'LineStyle','none','facealpha',.3);


hold on


plot(CircleVector,CrossCorrelation_mean,'color',[1 0 0],'LineWidth',1);
hold on

x2 = [CircleVector, fliplr(CircleVector)];
inBetween = [CrossCorrelation_mean-CrossCorrelation_std, fliplr(CrossCorrelation_mean+CrossCorrelation_std)];

fill(x2, inBetween, [1 0 0],'LineStyle','none','facealpha',.3);
yticks([-0.25:0.25:1.25]);
ylim([-0.25 1.25]);
xlim([0 18]);
xticks([0:2:18]);
box off
set(gca, 'TickDir', 'out')
