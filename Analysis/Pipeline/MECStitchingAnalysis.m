close all
clear all
clc

%% load in the color map
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
%%
SNR_threshold=3;
EventCount_threshold=100;
Overlapping_threshold=0.75;
%%
FOV_Total=5;
StichingPoor=cell(FOV_Total,1);
for    i=1:1:FOV_Total
    disp(['please select DataSet for FOV ', num2str(i)]);
    File_folder_raw=uigetdir;
    NATFilePath=[File_folder_raw,'\NAT.mat'];
    InformationFilePath=[File_folder_raw,'\ExperimentInformation.mat'];
    GridCellAnalysisFilePath=[File_folder_raw,'\GridCellAnalysis.mat'];
    StichingPoor{i,1}.NAT=load(NATFilePath);
    StichingPoor{i,1}.Information=load(InformationFilePath);
    StichingPoor{i,1}.GridCellAnalysis=load(GridCellAnalysisFilePath);
    disp(['loading DataSet for FOV ', num2str(i)]);
end
% %% reload the grid cells analysis
% for    i=1:1:FOV_Total
%     disp(['please select DataSet for FOV ', num2str(i)]);
%     File_folder_raw=uigetdir;
%     GridCellAnalysisFilePath=[File_folder_raw,'\GridCellAnalysis.mat'];
%     StichingPoor{i,1}.GridCellAnalysis=load(GridCellAnalysisFilePath);
%     disp(['loading DataSet for FOV ', num2str(i)]);
% end
%%
load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.mat');
FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.csv');

%% calculate total number of cells
TotalCell=0;
for i=1:1:FOV_Total
    TotalCell=StichingPoor{i,1}.Information.ExperimentInformation.TotalCell+TotalCell;
end
%%
FOVposition=load('FOVposition.mat');

% 
% FOVposition=[185,117;...
%     15,178;...
%     190,8;...
%     226,178;...
%     21,0];
%%
for    i=1:1:FOV_Total
    StichingPoor{i,1}.StitchingPosition=FOVposition(i,:);
end
%%
% FOV stitching and show channel 1
close all
figure
x0=10;
y0=10;
width=1400;
height=1400;
set(gcf,'position',[x0,y0,width,height])
%     a1=subplot(1,3,1,'align')
%     for FOV_ID=[2 1 3 4 5 ]
for FOV_ID= [2 1 4 5  3 ]
    
    Image_mean=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.ImagingOptions.meanImg(:,257:512);
    Image_mean_corrected=imwarp(Image_mean,TransformMatrix,'OutputView',imref2d(size(Image_mean)));
    X_shift=StichingPoor{FOV_ID,1}.StitchingPosition(1);
    Y_shift=StichingPoor{FOV_ID,1}.StitchingPosition(2);
    h=imagesc(X_shift,Y_shift,(Image_mean_corrected-min(min(Image_mean)))/(1*prctile(Image_mean(:),99)-min(min(Image_mean))),[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
    text(X_shift+128,Y_shift+128,['FOV',num2str(FOV_ID)],'Color','red','FontSize',14);
    hold on
    set(h,'AlphaData', 0.5*Image_mean_corrected>0);
    hold on
    colormap(gca,gray)
    daspect([1 1 1]);
    xlabel('x')
    ylabel('y')
    camroll(-90)
    %         set(gca,'XDir','reverse');
    %         set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
end
xlim([0 500]);
ylim([-10 450]);
set(gca,'Color',[0.5 0.5 0.5]);
set(gcf,'Color',[0.5 0.5 0.5]);
% grid on
set(gca, 'XColor', 'r')
xticks([0:2:500]);
yticks([0:2:500]);
%%
% FOV stitching and show channel 2
close all
figure
x0=10;
y0=10;
width=1400;
height=1400;
set(gcf,'position',[x0,y0,width,height])
%     a1=subplot(1,3,1,'align')
%     for FOV_ID=[2 1 3 4 5 ]
for FOV_ID= [2 1 4 5  3 ]
    
    Image_mean=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.ImagingOptions.meanImg_chan2(:,257:512);
    Image_mean_corrected=imwarp(Image_mean,TransformMatrix,'OutputView',imref2d(size(Image_mean)));
    X_shift=StichingPoor{FOV_ID,1}.StitchingPosition(1);
    Y_shift=StichingPoor{FOV_ID,1}.StitchingPosition(2);
    h=imagesc(X_shift,Y_shift,(Image_mean_corrected-min(min(Image_mean)))/(1*prctile(Image_mean(:),99)-min(min(Image_mean))),[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
    text(X_shift+128,Y_shift+128,['FOV',num2str(FOV_ID)],'Color','red','FontSize',14);
    hold on
    set(h,'AlphaData', 0.5*Image_mean_corrected>0);
    hold on
    colormap(gca,gray)
    daspect([1 1 1]);
    xlabel('x')
    ylabel('y')
    camroll(-90)
    %         set(gca,'XDir','reverse');
    %         set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
end
xlim([0 500]);
ylim([-10 450]);
set(gca,'Color',[0.5 0.5 0.5]);
set(gcf,'Color',[0.5 0.5 0.5]);
% grid on
set(gca, 'XColor', 'r')
xticks([0:2:500]);
yticks([0:2:500]);
%% plot all cells according to different FOV
close all
CMP=jet(5);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
k=1;
Plane=2;
for FOV_ID=1:1:5
    for i=1:1:StichingPoor{FOV_ID,1}.Information.ExperimentInformation.TotalCell
        if StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellinPlane(i)==Plane-1;
            Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
            Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
            Xpoint(Xpoint==0)=256;
            Ypoint(Ypoint==0)=256;
            ROI_pannel=zeros(size(Image_mean));
            for m=1:1:length(Xpoint)
                ROI_pannel(Ypoint(m),Xpoint(m))=1;
            end
            ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
            [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
            %                 Center_x(k+1)=mean(Xpoint)+FOVposition(FOV_ID,1);
            %                 Center_y(k+1)=mean(Ypoint)+FOVposition(FOV_ID,2);
            P=convhull(Xpoint,Ypoint);
            fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(FOV_ID,:),'LineStyle','none','facealpha',.6);
            hold on;
        else
        end
    end
    text(FOVposition(FOV_ID,1)+128,FOVposition(FOV_ID,2)+128,['FOV',num2str(FOV_ID)],'Color',CMP(FOV_ID,:),'FontSize',20);
end
hold off
box off
%     axis off
axis square
daspect([1 1 1]);
xlim([0 450]);
xlabel('x')
ylim([0 450]);
ylabel('y')
colormap(a3,gray)
caxis([0 1] );
camroll(-90)
set(gca,'YDir','reverse');
set(gcf,'color',[0 0 0]);
set(gca,'color',[0 0 0]);
set(a3, 'XAxisLocation', 'top')
set(a3, 'TickDir', 'out')
%% draw ROI of MEC beased on channel two image
ROI_MEC={};
disp('Now please drop the ROI for MEC')
ROI_MEC{1,1}=drawfreehand('Color',color_scheme_aaas(1,:));
%%
MEC_border=ROI_MEC{1,1}.Position;
%% find cells belong to MEC
% colume1: FOV ID;
% colume2: Cell ID;
% colume3: stitched x position
% colume4: stitched y position
% colume5: plane
% colume6: belong to MEC or not
% colume7: SNR
% colume8: event count
% colume9: Is grid cell or not
% colume10: grid score;
% colume11: Is repeared cell or not;
% colume12: Best event train shift;
% colume13: Filted grid cell;
% colume14~16: reserved;
NeuronMatrix=zeros(TotalCell,16);
k=1;
for FOV_ID=1:1:5
    for i=1:1:StichingPoor{FOV_ID,1}.Information.ExperimentInformation.TotalCell
        Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
        Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
        Xpoint(Xpoint==0)=256;
        Ypoint(Ypoint==0)=256;
        ROI_pannel=zeros(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.Framesize);
        for m=1:1:length(Xpoint)
            ROI_pannel(Ypoint(m),Xpoint(m))=1;
        end
        ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
        Center_x=mean(Xpoint);
        Center_y=mean(Ypoint);
        NeuronMatrix(k,1)=FOV_ID;
        NeuronMatrix(k,2)=i;
        NeuronMatrix(k,3)=FOVposition(FOV_ID,1)+Center_x;
        NeuronMatrix(k,4)=FOVposition(FOV_ID,2)+Center_y;
        NeuronMatrix(k,5)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellinPlane(i)+1;
        if inROI(ROI_MEC{1,1},NeuronMatrix(k,3),NeuronMatrix(k,4))==1
            NeuronMatrix(k,6)=1;
            
        else
            NeuronMatrix(k,6)=0;
        end
        NeuronMatrix(k,7)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellSNR(i);
        NeuronMatrix(k,8)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.EventCount_raw(i);
        if ismember(i,StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.IsGridCell{1,1})
            NeuronMatrix(k,9)=1;
        else
            NeuronMatrix(k,9)=0;
        end
        NeuronMatrix(k,10)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,end);
        if ismember(i,StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.Gridcell_overlapped)
            NeuronMatrix(k,11)=1;
        else
        end
        if ismember(i,StichingPoor{FOV_ID,1}.Information.ExperimentInformation.RepeatCell)
            NeuronMatrix(k,11)=1;
        else
        end
        NeuronMatrix(k,12)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridBestShift(i);
        k=k+1;
        disp(['Now checking cell ',num2str(k-1)]);
    end
end

%% filter out grid cell with SNR and events count
NeuronMatrix(:,13)=NeuronMatrix(:,9).*(NeuronMatrix(:,7)>SNR_threshold).*(NeuronMatrix(:,8)>EventCount_threshold).*(1-NeuronMatrix(:,11));
%% filter out  cell with overlapping
MaxDistance=10^2;
Cell_P1=find((NeuronMatrix(:,5)==1));
Cell_P2=find((NeuronMatrix(:,5)==2));
Overlap_P1=zeros(size(Cell_P1,1),size(Cell_P1,1),2);
Overlap_P2=zeros(size(Cell_P2,1),size(Cell_P2,1),2);
for k=1:1:size(Overlap_P1,2)
    FOV_ID=NeuronMatrix(Cell_P1(k),1);
    i=NeuronMatrix(Cell_P1(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.Framesize);
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    Xpoint=Xpoint+FOVposition(FOV_ID,1);
    Ypoint=Ypoint+FOVposition(FOV_ID,2);
    Center_x_1=mean(Xpoint);
    Center_y_1=mean(Ypoint);
    ROI_size_1=length(find(ROI_pannel_corrected(:)>0));
    for p=1:1:size(Overlap_P1,2)
        FOV_ID_2=NeuronMatrix(Cell_P1(p),1);
        i_2=NeuronMatrix(Cell_P1(p),2);
        Xpoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.xpix)',256);
        Ypoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.ypix)',256);
        Xpoint_2(Xpoint_2==0)=256;
        Ypoint_2(Ypoint_2==0)=256;
        ROI_pannel_2=zeros(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.Framesize);
        for m=1:1:length(Xpoint_2)
            ROI_pannel_2(Ypoint_2(m),Xpoint_2(m))=1;
        end
        ROI_pannel_2_corrected=imwarp(ROI_pannel_2,TransformMatrix,'OutputView',imref2d(size(ROI_pannel_2)));
        [Ypoint_2,Xpoint_2,~]=find(ROI_pannel_2_corrected>0);
        ROI_size_2=length(find(ROI_pannel_2_corrected(:)>0));
        Xpoint_2=Xpoint_2+FOVposition(FOV_ID_2,1);
        Ypoint_2=Ypoint_2+FOVposition(FOV_ID_2,2);
        Center_x_2=mean(Xpoint_2);
        Center_y_2=mean(Ypoint_2);
        D=(Center_x_2-Center_x_1)^2+(Center_y_2-Center_y_1)^2;
        if D<MaxDistance
            Overlappannel=zeros(600,600);
            for j=1:1:length(Xpoint)
                Overlappannel(Xpoint(j),Ypoint(j))=1;
            end
            for j=1:1:length(Xpoint_2)
                Overlappannel(Xpoint_2(j),Ypoint_2(j))=Overlappannel(Xpoint_2(j),Ypoint_2(j))+1;
            end
            Ratio_overlap1=length(find(Overlappannel(:)==2))./ROI_size_1;
            Ratio_overlap2=length(find(Overlappannel(:)==2))./ROI_size_2;
            Overlap_P1(k,p,1)=Ratio_overlap1;
            Overlap_P1(k,p,2)=Ratio_overlap2;
        else
        end
        disp([num2str(k),' vs ', num2str(p)]);
    end
end
%
for k=1:1:size(Overlap_P2,2)
    FOV_ID=NeuronMatrix(Cell_P2(k),1);
    i=NeuronMatrix(Cell_P2(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.Framesize);
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    Xpoint=Xpoint+FOVposition(FOV_ID,1);
    Ypoint=Ypoint+FOVposition(FOV_ID,2);
    ROI_size_1=length(find(ROI_pannel_corrected(:)>0));
    Center_x_1=mean(Xpoint)+FOVposition(FOV_ID,1);
    Center_y_1=mean(Ypoint)+FOVposition(FOV_ID,2);
    for p=1:1:size(Overlap_P2,2)
        Overlappannel=zeros(600,600);
        for j=1:1:length(Xpoint)
            Overlappannel(Xpoint(j),Ypoint(j))=1;
        end
        FOV_ID_2=NeuronMatrix(Cell_P2(p),1);
        i_2=NeuronMatrix(Cell_P2(p),2);
        Xpoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.xpix)',256);
        Ypoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.ypix)',256);
        Xpoint_2(Xpoint_2==0)=256;
        Ypoint_2(Ypoint_2==0)=256;
        ROI_pannel_2=zeros(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.Framesize);
        for m=1:1:length(Xpoint_2)
            ROI_pannel_2(Ypoint_2(m),Xpoint_2(m))=1;
        end
        ROI_pannel_2_corrected=imwarp(ROI_pannel_2,TransformMatrix,'OutputView',imref2d(size(ROI_pannel_2)));
        [Ypoint_2,Xpoint_2,~]=find(ROI_pannel_2_corrected>0);
        ROI_size_2=length(find(ROI_pannel_2_corrected(:)>0));
        Xpoint_2=Xpoint_2+FOVposition(FOV_ID_2,1);
        Ypoint_2=Ypoint_2+FOVposition(FOV_ID_2,2);
        Center_x_2=mean(Xpoint_2)+FOVposition(FOV_ID_2,1);
        Center_y_2=mean(Ypoint_2)+FOVposition(FOV_ID_2,2);
        D=(Center_x_2-Center_x_1)^2+(Center_y_2-Center_y_1)^2;
        if D<MaxDistance
            for j=1:1:length(Xpoint_2)
                Overlappannel(Xpoint_2(j),Ypoint_2(j))=Overlappannel(Xpoint_2(j),Ypoint_2(j))+1;
            end
            Ratio_overlap1=length(find(Overlappannel(:)==2))./ROI_size_1;
            Ratio_overlap2=length(find(Overlappannel(:)==2))./ROI_size_2;
            Overlap_P2(k,p,1)=Ratio_overlap1;
            Overlap_P2(k,p,2)=Ratio_overlap2;
        else
        end
        disp([num2str(k),' vs ', num2str(p)]);
    end
end
%%
for i=1:1: size(Overlap_P1,1)
    for j=1:1: size(Overlap_P1,2)
        MeanOverlap_P1(i,j)=max(Overlap_P1(i,j,:));
    end
end

for i=1:1: size(Overlap_P2,1)
    for j=1:1: size(Overlap_P2,2)
        MeanOverlap_P2(i,j)=max(Overlap_P2(i,j,:));
    end
end
%%
[RepeatCellPair_P1_C1 RepeatCellPair_P1_C2]=find((MeanOverlap_P1<1).*(MeanOverlap_P1>Overlapping_threshold)==1);
[RepeatCellPair_P2_C1 RepeatCellPair_P2_C2]=find((MeanOverlap_P2<1).*(MeanOverlap_P2>Overlapping_threshold)==1);
%% plot the overlapmap
close all
figure
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
imagesc(MeanOverlap_P2(:,:));
colormap(jet)
caxis([0 1] );
colorbar
daspect([1 1 1]);
%%
Tem=MeanOverlap_P1(:);
OverlapDistribution_P1=Tem(find((MeanOverlap_P1(:)<1).*(MeanOverlap_P1(:)>0)==1));
close all
figure
x0=10;
y0=10;
width=250;
height=350;
set(gcf,'position',[x0,y0,width,height])

histogram(OverlapDistribution_P1,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor',[1 1 1])


xlim([0 1])
xticks([0:0.25:1]);
box off
set(gca, 'TickDir', 'out')
%%
Tem=MeanOverlap_P2(:);
OverlapDistribution_P2=Tem(find((MeanOverlap_P2(:)<1).*(MeanOverlap_P2(:)>0)==1));
close all
figure
x0=10;
y0=10;
width=250;
height=350;
set(gcf,'position',[x0,y0,width,height])

histogram(OverlapDistribution_P2,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor',[1 1 1])


xlim([0 1])
xticks([0:0.25:1]);
box off
set(gca, 'TickDir', 'out')
%%
Remove=[];
for h=1:1:size(RepeatCellPair_P1_C1,1)
    SNR1=NeuronMatrix(Cell_P1(RepeatCellPair_P1_C1(h)),7);
    SNR2=NeuronMatrix(Cell_P1(RepeatCellPair_P1_C2(h)),7);
    [~,Remove(h)]=min([SNR1,SNR2]);
end
RepeatCell_P1=[];
for i=1:1:size(Remove,2)
    if Remove(i)==1
        RepeatCell_P1(i)=RepeatCellPair_P1_C1(i);
    else
        RepeatCell_P1(i)=RepeatCellPair_P1_C2(i);
    end
end



Remove=[];
for h=1:1:size(RepeatCellPair_P2_C1,1)
    SNR1=NeuronMatrix(Cell_P2(RepeatCellPair_P2_C1(h)),7);
    SNR2=NeuronMatrix(Cell_P2(RepeatCellPair_P2_C2(h)),7);
    [~,Remove(h)]=min([SNR1,SNR2]);
end
RepeatCell_P2=[];
for i=1:1:size(Remove,2)
    if Remove(i)==1
        RepeatCell_P2(i)=RepeatCellPair_P2_C1(i);
    else
        RepeatCell_P2(i)=RepeatCellPair_P2_C2(i);
    end
end
RepeatCell_P1=unique(RepeatCell_P1,'first');
RepeatCell_P2=unique(RepeatCell_P2,'first');
%% plot all cells according to different FOV (remove repeated cells)
close all
CMP=jet(5);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
k=1;
Plane=1;
for h=1:1:size(Cell_P1,1)
    FOV_ID=NeuronMatrix(Cell_P1(h),1);
    i=NeuronMatrix(Cell_P1(h),2);
    if ~ismember(h,RepeatCell_P1)
        Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
        Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
        Xpoint(Xpoint==0)=256;
        Ypoint(Ypoint==0)=256;
        ROI_pannel=zeros(size(Image_mean));
        for m=1:1:length(Xpoint)
            ROI_pannel(Ypoint(m),Xpoint(m))=1;
        end
        ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
        P=convhull(Xpoint,Ypoint);
        fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(FOV_ID,:),'LineStyle','none','facealpha',.6);
        hold on;
    else
    end
    %         text(FOVposition(FOV_ID,1)+128,FOVposition(FOV_ID,2)+128,['FOV',num2str(FOV_ID)],'Color',CMP(FOV_ID,:),'FontSize',20);
end
hold off
box off
%     axis off
axis square
daspect([1 1 1]);
xlim([0 450]);
xlabel('x')
ylim([0 450]);
ylabel('y')
colormap(a3,gray)
caxis([0 1] );
camroll(-90)
set(gca,'YDir','reverse');
set(gcf,'color',[0 0 0]);
set(gca,'color',[0 0 0]);
set(a3, 'XAxisLocation', 'top')
set(a3, 'TickDir', 'out')

%%
CellinMEC=find(NeuronMatrix(:,6)==1);
CellinP1=find(NeuronMatrix(:,5)==1);
CellinP2=find(NeuronMatrix(:,5)==2);
CellinMECinP1=CellinMEC(ismember(CellinMEC,CellinP1));
CellinMECinP2=CellinMEC(ismember(CellinMEC,CellinP2));
CellinP1_noRepeat=CellinP1;
CellinP1_noRepeat(RepeatCell_P1)=[];
CellinP2_noRepeat=CellinP2;
CellinP2_noRepeat(RepeatCell_P2)=[];
CellinMEC_P1_nonrepeat=CellinMEC(ismember(CellinMEC,CellinP1_noRepeat));
CellinMEC_P2_nonrepeat=CellinMEC(ismember(CellinMEC,CellinP2_noRepeat));
%%
GridCell_total=find(NeuronMatrix(:,9)==1);
GridCell_filtered=find(NeuronMatrix(:,13)==1);
GridCell_P1_noRepeat=GridCell_filtered(ismember(GridCell_filtered,CellinP1_noRepeat));
GridCell_P2_noRepeat=GridCell_filtered(ismember(GridCell_filtered,CellinP2_noRepeat));
GridCell_P1_MEC_noRepeat=GridCell_P1_noRepeat(ismember(GridCell_P1_noRepeat,CellinMEC_P1_nonrepeat));
GridCell_P2_MEC_noRepeat=GridCell_P2_noRepeat(ismember(GridCell_P2_noRepeat,CellinMEC_P2_nonrepeat));
GridCell_inTotal_nonrepeat=[GridCell_P1_noRepeat;GridCell_P2_noRepeat];
GridCell_inMEC=[GridCell_P1_MEC_noRepeat;GridCell_P2_MEC_noRepeat];
%% show all grid cells according to their DV position
%
GridCellMatrix=NeuronMatrix(GridCell_inTotal_nonrepeat,:);

[~,DVrank_GC]=sort(GridCellMatrix(:,10),'descend');
%%
% GridCellMatrix=NeuronMatrix(GridCell_filtered,:);
%
% [~,DVrank_GC]=sort(GridCellMatrix(:,3),'ascend');

%% show spatial tuning maps for all grid cells
close all
Shuffling=1000;
% CellPerPlot=50;
k=1;
% Xrange=[0 500];
session=1;
figure
x0=10;
y0=10;
width=1000;
height=2000;
set(gcf,'position',[x0,y0,width,height])
j=1;
for k=1:1:size(GridCellMatrix,1)
    subplot(21,15,k,'align')
    FOV_ID=GridCellMatrix(DVrank_GC(k),1);
    i=GridCellMatrix(DVrank_GC(k),2);
    MAP=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.ActivityMap{1,i}.z;
    MAX=prctile(MAP(:),99.5);
    %             imagesc(flipud((MAP-min(min(MAP)))./(max(max(MAP))-min(min(MAP)))),'AlphaData',MAP>0);
    imagesc(flipud(MAP./MAX),'AlphaData',(flipud(MAP)>0));
    CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    %             colormap(parula)
    %             colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP,1)])
    xlim([0 size(MAP,2)])
    %             title(['#',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    %             title(['#',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    %             title(['P ',num2str(max(MAP(:)),'%.2f')],'FontSize',7);
    
    daspect([1 1 1]);
    box off
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);
%% show autocorrelation maps for all grid cells
close all
Shuffling=1000;
% CellPerPlot=50;
k=1;
% Xrange=[0 500];
session=1;
figure
x0=10;
y0=10;
width=1000;
height=2000;
set(gcf,'position',[x0,y0,width,height])
j=1;
for k=1:1:size(GridCellMatrix,1)
    subplot(21,15,k,'align')
    FOV_ID=GridCellMatrix(DVrank_GC(k),1);
    i=GridCellMatrix(DVrank_GC(k),2);
    MAP=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.AutocorrelationMap{1,i};
    MIN=min(MAP(:));
    %             imagesc(flipud((MAP-min(min(MAP)))./(max(max(MAP))-min(min(MAP)))),'AlphaData',MAP>0);
    imagesc(flipud((MAP-MIN)./(1-MIN)));
    %             CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    %             colormap(parula)
    colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP,1)])
    xlim([0 size(MAP,2)])
    %             title(['#',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    %             title(['G ',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f')]);
    
    daspect([1 1 1]);
    box off
    axis off
end
set(gca,'color',[1 1 1]);
set(gcf,'color',[1 1 1]);

%% show events plot for all grid cells
close all
Shuffling=1000;
% CellPerPlot=50;
k=1;
% Xrange=[0 500];
session=1;
figure
x0=-100;
y0=-100;
width=1000;
height=2000;
set(gcf,'position',[x0,y0,width,height])
j=1;
coluum=15;
Shift=42;
S=0;
CMP=winter(50);
for k=1:1:size(GridCellMatrix,1)
    P=mod(k,coluum);
    if P==0
        P=coluum;
        
    elseif P==1
        S=S+1;
    else
    end
    
    %         subplot(20,15,k,'align')
    FOV_ID=GridCellMatrix(DVrank_GC(k),1);
    i=GridCellMatrix(DVrank_GC(k),2);
    Bestshift=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridBestShift(i);
    SelectedFrame_filtered=intersect(find(~isnan(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+10))),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    Event_filtered=intersect(find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+12)>0),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
    Event_filtered=intersect(Event_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
    
    Position=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(Event_filtered,2:3);
    Position(:,1)=circshift(Position(:,1),Bestshift);
    Position(:,2)=circshift(Position(:,2),Bestshift);
    Event=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(Event_filtered,4*i+12);
    Max=max(Event(:));
    %             scatter(Position(:,1)+85*(P-1),Position(:,2)-85*(S-1),8*Event./Max,CMP(51-S,:),'filled',0.5)
    %             scatter(Position(:,1)+90*(P-1),Position(:,2)-90*(S-1),3*sqrt(sqrt(Event./Max)),CMP(51-S,:),'filled','MarkerFaceAlpha',0.8)
    scatter(Position(:,1)+92*(P-1),Position(:,2)-92*(S-1),3*sqrt(sqrt(Event./Max)),[1 0 0],'filled','MarkerFaceAlpha',0.8)
    %             text(92*(P-1),-92*(S-1),[num2str(FOV_ID),'-',num2str(i)]);
    ylim([-92*21,48])
    xlim([-48,92*14])
    %             title(['#',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(max(max(MAP)),'%.2f')]);
    %             title(['#',num2str(i)]);
    
    daspect([1 1 1]);
    hold on
    
    axis off
end
box off
% set(gca,'color',[0 0 0]);
% set(gcf,'color',[0 0 0]);
%% plot thre SNR for all grid cells
figure
x0=100;
y0=100;
width=100;
height=300;
set(gcf,'position',[x0,y0,width,height])
SNR_gridcell=GridCellMatrix(:,7);
WJplots.Violinplot.violinplot(SNR_gridcell);

%%
GS_gridcell=zeros(size(GridCellMatrix,1),1003);

for k=1:1:size(GridCellMatrix,1)
    FOV_ID=GridCellMatrix(DVrank_GC(k),1);
    i=GridCellMatrix(DVrank_GC(k),2);
    GS_gridcell(k,:)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,:);
end
%%
GS_shuffled=GS_gridcell(:,1:1000);
close all
figure
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])

yyaxis left
histogram(GS_shuffled(:),'Normalization','probability','BinWidth',0.05,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[1 1 1]);
set(gca,  'XColor','k', 'YColor','k')
ylim([0 0.15])
yyaxis right
histogram(GS_gridcell(:,1003),'Normalization','probability','BinWidth',0.1,'DisplayStyle','stairs','LineWidth',2);
ylim([0 0.15])
xlim([-1,1.8])
xticks([-1:0.5:2])
set(gca, 'TickDir', 'out')
box off
%% chance level analysis
GSinMEC_shuffled=zeros(size(GridCell_inMEC,1),1003);
for k=1:1:size(GridCell_inMEC,1)
    FOV_ID=NeuronMatrix(GridCell_inMEC(k),1);
    i=NeuronMatrix(GridCell_inMEC(k),2);
    GSinMEC_shuffled(k,:)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,:);
end

%% Real GCs v.s. change
j=1; %session ID
TotalRun=200;
TotalNumber=length(CellinMEC);
ChanceWin=zeros(TotalRun,1);
block=ceil(TotalNumber./size(GridCell_inMEC,1));
Blockvector=[1:block:block*(TotalRun+1)];
%%
for w=1:1:TotalRun
    TotalNumber_GC_chance=0;
    while (TotalNumber_GC_chance<TotalNumber)
        for m=1:1:block
            Shullfed=GSinMEC_shuffled(:,[Blockvector(w)+m-1,1001]);
            for i=1:1:size(Shullfed,1)
                if Shullfed(i,1)>Shullfed(i,2)
                    ChanceWin(w,1)= ChanceWin(w,1)+1;
                else
                end
            end
        end
        TotalNumber_GC_chance=TotalNumber_GC_chance+block*size(GSinMEC_shuffled,1);
    end
end
%%
close all
figure
x0=100;
y0=100;
width=280;
height=330;
set(gcf,'position',[x0,y0,width,height])
histogram(100*ChanceWin/TotalNumber_GC_chance,10,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
hold
plot(100.*[size(GridCell_inMEC,1)/TotalNumber size(GridCell_inMEC,1)/TotalNumber],[0 0.15],'LineWidth',3,'color',[1 0 0]);
box off
set(gca, 'TickDir', 'out')
xlim([-0.5 40])
xticks([0:5:40])
ylim([0 0.25])
%% plot all cells with grid cells
close all
CMP=jet(5);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
for k=1:1:size(CellinP1_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP1_noRepeat(k),1);
    i=NeuronMatrix(CellinP1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[.3 .3 .3],'LineStyle','none','facealpha',.3);
    hold on;
end

for k=1:1:size(CellinP2_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP2_noRepeat(k),1);
    i=NeuronMatrix(CellinP2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[.3 .3 .3],'LineStyle','none','facealpha',.3);
    hold on;
end

for k=1:1:size(GridCell_P2_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P2_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[102 153 255]/255,'LineStyle','none','facealpha',.9);
    hold on;
end


for k=1:1:size(GridCell_P1_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P1_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),color_scheme_aaas(2,:),'LineStyle','none','facealpha',.9);
    hold on;
end

hold off
box off
%     axis off
axis square
daspect([1 1 1]);
xlim([0 450]);
xlabel('x')
ylim([0 450]);
ylabel('y')
%     colormap(a3,gray)
caxis([0 1] );
camroll(-90)
set(gca,'YDir','reverse');
set(gcf,'color',[0 0 0]);
set(gca,'color',[0 0 0]);
set(a3, 'XAxisLocation', 'top')
set(a3, 'TickDir', 'out')

axis off
box off
%% plot a example Grid Cell
close all
CMP=WJplots.CMP.inferno(256);
Shuffling=1000;
figure
width=360;
height=360;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
j=1;
FOV_ID=1;
i=46;
Peak=max(max(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.ActivityMap{1,i}.z));
subplot(2,2,1)
Bestshift=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridBestShift(i);
SelectedFrame_filtered=intersect(find(~isnan(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+10))),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
Event_filtered=intersect(find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+12)>0),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
Event_filtered=intersect(Event_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>2.5));% filter out the frames with speed threadhold
AnimaPosition=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(SelectedFrame_filtered,2:3);
EventPosition=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(Event_filtered,2:3);
EventPosition(:,1)=circshift(EventPosition(:,1),Bestshift);
EventPosition(:,2)=circshift(EventPosition(:,2),Bestshift);
Event=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(Event_filtered,4*i+12);
Max=max(Event(:));
%             scatter(Position(:,1)+85*(P-1),Position(:,2)-85*(S-1),8*Event./Max,CMP(51-S,:),'filled',0.5)
%             scatter(Position(:,1)+90*(P-1),Position(:,2)-90*(S-1),3*sqrt(sqrt(Event./Max)),CMP(51-S,:),'filled','MarkerFaceAlpha',0.8)
plot(AnimaPosition(:,1),AnimaPosition(:,2),'color', [0.6 0.6 0.6],'LineWidth',0.5);
hold on
scatter(EventPosition(:,1),EventPosition(:,2),40*(Event./Max),CMP(120,:),'filled','MarkerFaceAlpha',0.8)
ylim([-41,41])
xlim([-41,41])
title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
daspect([1 1 1]);
box off
axis off

subplot(2,2,2)
MAP=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.ActivityMap{1,i}.z;
MAX=prctile(MAP(:),99.5);
imagesc(flipud(MAP./MAX),'AlphaData',MAP>0);

colormap(CMP)
caxis([0 1] );
ylim([0 size(MAP,1)])
xlim([0 size(MAP,2)])
title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
daspect([1 1 1]);
box off
axis off


subplot(2,2,3)

% AutocorrelationMap=analyses.autocorrelation(MAP);
AutoMap=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.AutocorrelationMap{1,i};
imagesc(flipud((AutoMap+1)./2));
colormap(gca,jet)
caxis([0 1] );
ylim([0 size(AutoMap,1)])
xlim([0 size(AutoMap,2)])
title([num2str(FOV_ID),'-',num2str(i),' GC:',num2str(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,Shuffling+3,1),'%.2f'),' P:',num2str(Peak,'%.2f')]);
daspect([1 1 1]);
box off
axis off


subplot(2,2,4)

% AutocorrelationMap=analyses.autocorrelation(MAP);
Shuffled_example=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,1:1000);
H=histogram(Shuffled_example,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
hold on
GC=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,1003);
Cut=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,1001);
plot([GC GC],[0 max(H.Values)],'color',[1 0 0],'LineWidth',3);
hold on
plot([Cut Cut],[0 max(H.Values)],'color',[0 0.7 0.7],'LineWidth',3);
box off
ylim([0 ceil(max(H.Values)/0.05)*0.05])
set(gca, 'TickDir', 'out')


%% plot all cells with grid cells color coded by spacing
close all
CMP=jet(256);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
for k=1:1:size(CellinP1_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP1_noRepeat(k),1);
    i=NeuronMatrix(CellinP1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[0.3 0.3 0.3],'LineStyle','none','facealpha',.2);
    hold on;
end

for k=1:1:size(CellinP2_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP2_noRepeat(k),1);
    i=NeuronMatrix(CellinP2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[0.3 0.3 0.3],'LineStyle','none','facealpha',.3);
    hold on;
end

for k=1:1:size(GridCell_P2_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P2_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    Spacing=mean(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.spacing);
    Spacing=Spacing*2;
    Color=round((Spacing-20)/(60-20)*256);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(Color,:),'LineStyle','none','facealpha',.9);
    hold on;
end


for k=1:1:size(GridCell_P1_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P1_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    Spacing=mean(StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.spacing);
    Spacing=Spacing*2;
    Color=round((Spacing-20)/(60-20)*256);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(Color,:),'LineStyle','none','facealpha',.9);
    hold on;
end
hold off
box off
%     axis off
axis square
daspect([1 1 1]);
xlim([0 450]);
xlabel('x')
ylim([0 450]);
ylabel('y')
xticks([])
yticks([])
%     colormap(a3,gray)
caxis([0 1] );
camroll(-90)
set(gcf, 'InvertHardcopy', 'off')
set(gca,'YDir','reverse');
set(gcf,'color',[0 0 0]);
set(gca,'color',[0 0 0]);
set(a3, 'XAxisLocation', 'top')
set(a3, 'TickDir', 'out')
%     axis off
%     box off
%% calculate the grid cell offset
for k=1:1:size(GridCell_P1_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P1_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P1_noRepeat(k),2);
    Orientaition=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.orientation;
    Orientaition=[Orientaition;Orientaition+180];
    GC_P1_Offset(k)=preprocessing.DetectOffset(Orientaition,90);
end
for k=1:1:size(GridCell_P2_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P2_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P2_noRepeat(k),2);
    Orientaition=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.orientation;
    Orientaition=[Orientaition;Orientaition+180];
    GC_P2_Offset(k)=preprocessing.DetectOffset(Orientaition,90);
end

%% plot all cells with grid cells color coded by orientation
close all
CMP=jet(256);
figure(2)
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
for k=1:1:size(CellinP1_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP1_noRepeat(k),1);
    i=NeuronMatrix(CellinP1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[0.3 0.3 0.3],'LineStyle','none','facealpha',.2);
    hold on;
end

for k=1:1:size(CellinP2_noRepeat,1)
    FOV_ID=NeuronMatrix(CellinP2_noRepeat(k),1);
    i=NeuronMatrix(CellinP2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),[0.3 0.3 0.3],'LineStyle','none','facealpha',.3);
    hold on;
end

for k=1:1:size(GridCell_P2_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P2_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P2_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    Orientaition=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.orientation;
    Orientaition=[Orientaition;Orientaition+180];
    Offset=preprocessing.DetectOffset(Orientaition,90);
    if abs(Offset)<45
        Color=round((Offset+45)/(90)*255);
    else
        Color=128;
    end
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(Color,:),'LineStyle','none','facealpha',.9);
    hold on;
end


for k=1:1:size(GridCell_P1_noRepeat,1)
    FOV_ID=NeuronMatrix(GridCell_P1_noRepeat(k),1);
    i=NeuronMatrix(GridCell_P1_noRepeat(k),2);
    Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
    Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
    Xpoint(Xpoint==0)=256;
    Ypoint(Ypoint==0)=256;
    ROI_pannel=zeros(size(Image_mean));
    for m=1:1:length(Xpoint)
        ROI_pannel(Ypoint(m),Xpoint(m))=1;
    end
    ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));
    [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);
    P=convhull(Xpoint,Ypoint);
    Orientaition=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridsStat{1,i}.orientation;
    Orientaition=[Orientaition;Orientaition+180];
    Offset=preprocessing.DetectOffset(Orientaition,90);
    if abs(Offset)<45
        Color=round((Offset+45)/(90)*255);
    else
        Color=128;
    end
    fill(FOVposition(FOV_ID,1)+Xpoint(P),FOVposition(FOV_ID,2)+Ypoint(P),CMP(Color,:),'LineStyle','none','facealpha',.9);
    hold on;
end
hold off
box off
%     axis off
axis square
daspect([1 1 1]);
xlim([0 450]);
xlabel('x')
ylim([0 450]);
ylabel('y')
%     colormap(a3,gray)
caxis([0 1] );
camroll(-90)
set(gca,'YDir','reverse');
set(gcf,'color',[0 0 0]);
set(gca,'color',[0 0 0]);
set(a3, 'XAxisLocation', 'top')
set(a3, 'TickDir', 'out')
xticks([])
yticks([])
%     axis off
%     box off

box off
