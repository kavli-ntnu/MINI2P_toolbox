close all
clear all
clc
%%
SNR_threshold=3;
EventCount_threshold=100;
Overlapping_threshold=0.75;
%%
FOV_Total=2;
StichingPoor=cell(FOV_Total,1);
for    i=1:1:FOV_Total
    disp(['please select DataSet for FOV ', num2str(i)]);
    File_folder_raw=uigetdir;
    NAAKFilePath=[File_folder_raw,'\NAAK.mat'];  
    InformationFilePath=[File_folder_raw,'\ExperimentInformation.mat'];
    StichingPoor{i,1}.NAAK=load(NAAKFilePath);
    StichingPoor{i,1}.Information=load(InformationFilePath);
    disp(['loading DataSet for FOV ', num2str(i)]);
end
%%
load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.mat');
FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\OPENMINI2P2_Workshop\LFOV001\AVG_OpenMINI2P_LFOV001_50umGrind_256_00001_2021-03-10-16-16-17.csv');                     

%% calculate total number of cells
TotalCell=0;
    for i=1:1:FOV_Total
        TotalCell=StichingPoor{i,1}.Information.ExperimentInformation.TotalCell+TotalCell;
    end
%%   
FOVposition=[10,0;...
             0,163 ]+50;
%%
for    i=1:1:FOV_Total
    StichingPoor{i,1}.StitchingPosition=FOVposition(i,:);
end 
% FOV to check
%
FOVtoCheck=[1 2];
%
% step1:primary FOV alignment by overlapping landmarks
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
    for FOV_ID=FOVtoCheck     
        Image_mean=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.ImagingOptions.meanImg(:,257:512);
        Image_mean_corrected=imwarp(Image_mean,TransformMatrix,'OutputView',imref2d(size(Image_mean)));  
        Image_mean_corrected= imrotate(Image_mean_corrected,90);
        Image_mean_corrected = flip(Image_mean_corrected,2);
        X_shift=StichingPoor{FOV_ID,1}.StitchingPosition(1);
        Y_shift=StichingPoor{FOV_ID,1}.StitchingPosition(2); 
        h=imagesc(X_shift,Y_shift,(Image_mean_corrected-min(min(Image_mean)))/(1*prctile(Image_mean(:),99)-min(min(Image_mean))),'AlphaData', 0.1*Image_mean_corrected>0);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
        text(X_shift+128,Y_shift+128,['FOV',num2str(FOV_ID)],'Color','red','FontSize',14);
        hold on
%         set(h,'AlphaData', 0.1*Image_mean_corrected>0);
%         hold on
        colormap(gca,gray)
        daspect([1 1 1]);
        xlabel('x')
        ylabel('y')

%         set(gca,'XDir','reverse');
%         set(gca,'YDir','reverse');
        set(gca, 'TickDir', 'out')   
    end
    
% camroll(90)
% set(gca,'YDir','reverse');
% set(gca,'XDir','reverse');
xlim([0 500]);
ylim([0 500]);
set(gca,'Color',[0.5 0.5 0.5]);
set(gcf,'Color',[0.5 0.5 0.5]);
% grid on
set(gca, 'XColor', 'r')
xticks([0:2:500]);
yticks([0:2:500]);
%% (option) FOV stitching and show channel 2
% % 
% close all
% figure
% x0=10;
% y0=10;
% width=1400;
% height=1400;
% set(gcf,'position',[x0,y0,width,height])
% %     a1=subplot(1,3,1,'align')
% %     for FOV_ID=[2 1 3 4 5 ]
%     for FOV_ID= [2 1 4 5  3 ]     
% 
%         Image_mean=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.ImagingOptions.meanImg_chan2(:,257:512);
%         Image_mean_corrected=imwarp(Image_mean,TransformMatrix,'OutputView',imref2d(size(Image_mean)));  
%         X_shift=StichingPoor{FOV_ID,1}.StitchingPosition(1);
%         Y_shift=StichingPoor{FOV_ID,1}.StitchingPosition(2);
%         h=imagesc(X_shift,Y_shift,(Image_mean_corrected-min(min(Image_mean)))/(1*prctile(Image_mean(:),99)-min(min(Image_mean))),[0 1]);% ,'AlphaData',FullNeuronBehaviorDataSet.CellMasks{:,:,i}
%         text(X_shift+128,Y_shift+128,['FOV',num2str(FOV_ID)],'Color','red','FontSize',14);
%         hold on
%         set(h,'AlphaData', 0.5*Image_mean_corrected>0);
%         hold on
%         colormap(gca,gray)
%         daspect([1 1 1]);
%         xlabel('x')
%         ylabel('y')
%         camroll(-90)
% %         set(gca,'XDir','reverse');
% %         set(gca,'YDir','reverse');
%         set(gca, 'TickDir', 'out')   
%     end
% xlim([0 500]);
% ylim([-10 450]);
% set(gca,'Color',[0.5 0.5 0.5]);
% set(gcf,'Color',[0.5 0.5 0.5]);
% % grid on
% set(gca, 'XColor', 'r')
% xticks([0:2:500]);
% yticks([0:2:500]);
%% Step2: Refinement of the alignment by monitoring the image correlation
close all
Scanrange_x=[-3:1:7];
Scanrange_y=[-13:1:-3];
%
ImageCorrelation=zeros(length(Scanrange_x),length(Scanrange_y));
%
Image1_extend=nan(StichingPoor{FOVtoCheck(1),1}.Information.ExperimentInformation.Framesize*2);
Crop_x=1+FOVposition(1,1):StichingPoor{FOVtoCheck(1),1}.Information.ExperimentInformation.Framesize(1)+FOVposition(1,1);
Crop_y=1+FOVposition(1,2):StichingPoor{FOVtoCheck(1),1}.Information.ExperimentInformation.Framesize(2)+FOVposition(1,2);
ImageCorrect1=imwarp(StichingPoor{FOVtoCheck(1),1}.Information.ExperimentInformation.ImagingOptions.meanImg(:,257:512),TransformMatrix,'OutputView',imref2d(size(Image_mean)));  
ImageCorrect1= imrotate(ImageCorrect1,90);
ImageCorrect1 = flip(ImageCorrect1,2);
Image1_extend(Crop_y,Crop_x)=ImageCorrect1;% check the correlation of P1

%
% hold on
% imagesc(Image1_extend)
% daspect([1 1 1]);
% %
n=1;
for i=1:1:size(ImageCorrelation,1)   
    for j=1:1:size(ImageCorrelation,2)
        Image2_extend=nan(size(Image1_extend));
        Crop_x_2=1+FOVposition(2,1)+Scanrange_x(i):StichingPoor{FOVtoCheck(2),1}.Information.ExperimentInformation.Framesize(1)+FOVposition(2,1)+Scanrange_x(i);
        Crop_y_2=1+FOVposition(2,2)+Scanrange_y(j):StichingPoor{FOVtoCheck(2),1}.Information.ExperimentInformation.Framesize(2)+FOVposition(2,2)+Scanrange_y(j);
        ImageCorrect2=imwarp(StichingPoor{FOVtoCheck(2),1}.Information.ExperimentInformation.ImagingOptions.meanImg(:,257:512),TransformMatrix,'OutputView',imref2d(size(Image_mean)));  
        ImageCorrect2= imrotate(ImageCorrect2,90);
        ImageCorrect2 = flip(ImageCorrect2,2);                
        Image2_extend(Crop_y_2,Crop_x_2)=ImageCorrect2;       
        Overlap=zeros(size(Image1_extend));
            for p=1:1:size(Image1_extend,1)
                for k=1:1:size(Image1_extend,2)
                    if ~isnan(Image1_extend(p,k)) && ~isnan(Image2_extend(p,k))
                        Overlap(p,k)=1;
                    else    
                    end
                end
            end
        [X Y]=find(Overlap==1);
        Min_x=min(X);
        Max_x=max(X);
        Min_y=min(Y);
        Max_y=max(Y);
        if Max_x>Min_x && Max_y>Min_y
            ImageCorrelation(i,j)=corr2(Image1_extend(Min_x:Max_x,Min_y:Max_y),Image2_extend(Min_x:Max_x,Min_y:Max_y));     
        else
            ImageCorrelation(i,j)=0;
        end
%         Image3_extend(Crop_y,Crop_x)=ImageCorrect1;
%         Image3_extend(Crop_y_2,Crop_x_2)=ImageCorrect1;
%         subplot(6,ceil(length(CorrelationGrid(:))/6),k)
%         imagesc(Image3_extend)
        n=n+1
    end
end
%%
[Shift_x_P1, Shift_y_P2]=find(ImageCorrelation==max(ImageCorrelation(:)));
Best_shiftX=Scanrange_x(Shift_x_P1);
Best_shiftY=Scanrange_x(Shift_y_P2);
FOVposition(2,1)=FOVposition(2,1)+Best_shiftX;
FOVposition(2,2)=FOVposition(2,2)+Best_shiftY;
% close all
%%
ImageCorrelation_correct=0.95*(ImageCorrelation-min(ImageCorrelation(:)))./(max(ImageCorrelation(:))-min(ImageCorrelation(:)));
close all
figure
imagesc(ImageCorrelation_correct)
 daspect([1 1 1]);
xlabel('x')
ylabel('y')    
xticks([1 6 11])  
yticks([1 6 11]) 
set(gca, 'TickDir', 'out')
caxis([0 1])
% colormap(jet)
colorbar
% ylim([5,40]);
% xlim([35,60]);

%% Final check of the alignment by repeat cell overlapping
% plot all cells according to different FOV
%refine stiching
FOVposition_correct(1,1)=FOVposition(1,1);
FOVposition_correct(1,2)=FOVposition(1,2);
FOVposition_correct(2,1)=FOVposition(2,1)-6;
FOVposition_correct(2,2)=FOVposition(2,2)-2;
%
close all
CMP=jet(5);
figure
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
k=1;
Plane=2;
CMP=[0,1,0;1,0,1];
    for FOV_ID=FOVtoCheck
% for FOV_ID=2
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
                ROI_pannel_corrected= imrotate(ROI_pannel_corrected,90);
                ROI_pannel_corrected = flip(ROI_pannel_corrected,2);
                if FOV_ID==2
                ROI_pannel_corrected = imrotate(ROI_pannel_corrected,0.6); 
                else
                end
                [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);      
%                 Center_x(k+1)=mean(Xpoint)+FOVposition(FOV_ID,1);
%                 Center_y(k+1)=mean(Ypoint)+FOVposition(FOV_ID,2);
                P=convhull(Xpoint,Ypoint);             
                fill(FOVposition_correct(FOV_ID,1)+Xpoint(P),FOVposition_correct(FOV_ID,2)+Ypoint(P),CMP(FOV_ID,:),'LineStyle','none','facealpha',.6);                    
                hold on;    
            else
            end
        end
%         text(FOVposition_correct(FOV_ID,1)+128,FOVposition_correct(FOV_ID,2)+128,['FOV',num2str(FOV_ID)],'Color',CMP(FOV_ID,:),'FontSize',20);
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
%     camroll(-90)
    set(gca,'YDir','reverse');
    set(gcf,'color',[0 0 0]);
    set(gca,'color',[0 0 0]);
    set(a3, 'XAxisLocation', 'top')
    set(a3, 'TickDir', 'out')
box off
axis off

%% Build neuronmitrix
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
for FOV_ID=1:1:2
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
            ROI_pannel_corrected= imrotate(ROI_pannel_corrected,90);
            ROI_pannel_corrected = flip(ROI_pannel_corrected,2);
            if FOV_ID==2
                ROI_pannel_corrected = imrotate(ROI_pannel_corrected,0.6); 
            else
            end
            [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);      
            Center_x=mean(Xpoint);
            Center_y=mean(Ypoint); 
%     
          
            NeuronMatrix(k,1)=FOV_ID;
            NeuronMatrix(k,2)=i;
            NeuronMatrix(k,3)=FOVposition_correct(FOV_ID,1)+Center_x;
            NeuronMatrix(k,4)=FOVposition_correct(FOV_ID,2)+Center_y;
            NeuronMatrix(k,5)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellinPlane(i)+1;
%             if inROI(ROI_MEC{1,1},NeuronMatrix(k,3),NeuronMatrix(k,4))==1
%                 NeuronMatrix(k,6)=1;
%   
%             else
%                 NeuronMatrix(k,6)=0;
%             end
            NeuronMatrix(k,7)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellSNR(i);
            NeuronMatrix(k,8)=StichingPoor{FOV_ID,1}.Information.ExperimentInformation.EventCount_raw(i);
%             if ismember(i,StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.IsGridCell{1,1})
%                 NeuronMatrix(k,9)=1;
%             else
%                 NeuronMatrix(k,9)=0;
%             end
%             NeuronMatrix(k,10)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridScore_shuffled(i,end);
%             if ismember(i,StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.Gridcell_overlapped)
%                 NeuronMatrix(k,11)=1;
%             else 
%             end
%             if ismember(i,StichingPoor{FOV_ID,1}.Information.ExperimentInformation.RepeatCell)
%                 NeuronMatrix(k,11)=1;
%             else
%             end 
%             NeuronMatrix(k,12)=StichingPoor{FOV_ID,1}.GridCellAnalysis.GridCellAnalysis.GridBestShift(i);
            k=k+1;
            disp(['Now checking cell ',num2str(k-1)]);
        end
end

%% filter out  cell with overlapping
    MaxDistance=10^2;
    Cell_P1=find((NeuronMatrix(:,5)==1));
    Cell_P2=find((NeuronMatrix(:,5)==2)); 
    Overlap_P1=zeros(size(Cell_P1,1),size(Cell_P1,1),2);
    Overlap_P2=zeros(size(Cell_P2,1),size(Cell_P2,1),2);
%     for k=1:1:size(Overlap_P1,2)
%         FOV_ID=NeuronMatrix(Cell_P1(k),1);
%         i=NeuronMatrix(Cell_P1(k),2);
%         Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
%         Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
%         Xpoint(Xpoint==0)=256;
%         Ypoint(Ypoint==0)=256;
%         ROI_pannel=zeros(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.Framesize);
%         for m=1:1:length(Xpoint)
%             ROI_pannel(Ypoint(m),Xpoint(m))=1;
%         end
%         ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));             
%         [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0); 
%         Xpoint=Xpoint+FOVposition(FOV_ID,1);
%         Ypoint=Ypoint+FOVposition(FOV_ID,2);
%         Center_x_1=mean(Xpoint)+FOVposition(FOV_ID,1);
%         Center_y_1=mean(Ypoint)+FOVposition(FOV_ID,2); 
%   
%         ROI_size_1=length(find(ROI_pannel_corrected(:)>0));
%         
%         for p=1:1:size(Overlap_P1,2)
%             FOV_ID_2=NeuronMatrix(Cell_P1(p),1);
%             i_2=NeuronMatrix(Cell_P1(p),2);
%             Xpoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.xpix)',256);
%             Ypoint_2=mod(double(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.CellStat{1,i_2}.ypix)',256);
%             Xpoint_2(Xpoint_2==0)=256;
%             Ypoint_2(Ypoint_2==0)=256;
%             ROI_pannel_2=zeros(StichingPoor{FOV_ID_2,1}.Information.ExperimentInformation.Framesize);
%             for m=1:1:length(Xpoint_2)
%                 ROI_pannel_2(Ypoint_2(m),Xpoint_2(m))=1;
%             end
%             ROI_pannel_2_corrected=imwarp(ROI_pannel_2,TransformMatrix,'OutputView',imref2d(size(ROI_pannel_2)));             
%             [Ypoint_2,Xpoint_2,~]=find(ROI_pannel_2_corrected>0);
%             ROI_size_2=length(find(ROI_pannel_2_corrected(:)>0));
%             Xpoint_2=Xpoint_2+FOVposition(FOV_ID_2,1);
%             Ypoint_2=Ypoint_2+FOVposition(FOV_ID_2,2);
%             Center_x_2=mean(Xpoint_2)+FOVposition(FOV_ID_2,1);
%             Center_y_2=mean(Ypoint_2)+FOVposition(FOV_ID_2,2); 
%             D=(Center_x_2-Center_x_1)^2+(Center_y_2-Center_y_1)^2;
%             if D<MaxDistance
%                 Overlappannel=zeros(600,600);
%                 for j=1:1:length(Xpoint)
%                     Overlappannel(Xpoint(j),Ypoint(j))=1;
%                 end
%                 for j=1:1:length(Xpoint_2)
%                     Overlappannel(Xpoint_2(j),Ypoint_2(j))=Overlappannel(Xpoint_2(j),Ypoint_2(j))+1;
%                 end
%                 Ratio_overlap1=length(find(Overlappannel(:)==2))./ROI_size_1;
%                 Ratio_overlap2=length(find(Overlappannel(:)==2))./ROI_size_2;
%                 Overlap_P1(k,p,1)=Ratio_overlap1;
%                 Overlap_P1(k,p,2)=Ratio_overlap2;
%             else
%             end
%             disp([num2str(k),' vs ', num2str(p)]);
%         end
%     end
% %
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
        ROI_pannel_corrected= imrotate(ROI_pannel_corrected,90);
        ROI_pannel_corrected = flip(ROI_pannel_corrected,2);
        if FOV_ID==2
           ROI_pannel_corrected = imrotate(ROI_pannel_corrected,0.6); 
        else
        end
        [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0); 
        Xpoint=Xpoint+FOVposition_correct(FOV_ID,1);
        Ypoint=Ypoint+FOVposition_correct(FOV_ID,2);
        ROI_size_1=length(find(ROI_pannel_corrected(:)>0));
        Center_x_1=mean(Xpoint);
        Center_y_1=mean(Ypoint);            
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
            ROI_pannel_2_corrected= imrotate(ROI_pannel_2_corrected,90);
            ROI_pannel_2_corrected = flip(ROI_pannel_2_corrected,2);
                if FOV_ID_2==2
                   ROI_pannel_2_corrected = imrotate(ROI_pannel_2_corrected,0.6); 
                else
                end
            [Ypoint_2,Xpoint_2,~]=find(ROI_pannel_2_corrected>0);
            ROI_size_2=length(find(ROI_pannel_2_corrected(:)>0));
            Xpoint_2=Xpoint_2+FOVposition_correct(FOV_ID_2,1);
            Ypoint_2=Ypoint_2+FOVposition_correct(FOV_ID_2,2);
            Center_x_2=mean(Xpoint_2);
            Center_y_2=mean(Ypoint_2); 
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

%% plot the overlapmap
close all
    figure 
    x0=10;
    y0=10;
    width=500;
    height=500;
    set(gcf,'position',[x0,y0,width,height])  
imagesc((Overlap_P2(:,:,1)+Overlap_P2(:,:,1))/2);
colormap(jet)
caxis([0 1] ); 
colorbar
xticks
daspect([1 1 1]); 

%%

for i=1:1: size(Overlap_P2,1)
    for j=1:1: size(Overlap_P2,2)
        MeanOverlap(i,j)=max(Overlap_P2(i,j,:));
    end
end
%%

MeanOverlap_1=MeanOverlap(:);
OverlappingDistribution=MeanOverlap_1(find((MeanOverlap(:)<1).*(MeanOverlap(:)>0)==1));
[RepeatCellPair_1 RepeatCellPair_2]=find((MeanOverlap<1).*(MeanOverlap>0.75)==1);
close all
figure 
    x0=10;
    y0=10;
    width=400;
    height=350;
    set(gcf,'position',[x0,y0,width,height])

histogram(OverlappingDistribution,'Normalization', 'probability','FaceColor',[0.3 0.3 0.3],'EdgeColor',[1 1 1])


    xlim([0 1])
    xticks([0:0.25:1]);
    box off
    set(gca, 'TickDir', 'out')

%%
Remove=[];
for h=1:1:size(RepeatCellPair_1,1)
        SNR1=NeuronMatrix(Cell_P2(RepeatCellPair_1(h)),7);
        SNR2=NeuronMatrix(Cell_P2(RepeatCellPair_2(h)),7);
        [~,Remove(h)]=min([SNR1,SNR2]);
end
RepeatCell=[];
for i=1:1:size(Remove,2)
    if Remove(i)==1
        RepeatCell(i)=RepeatCellPair_1(i);
    else
        RepeatCell(i)=RepeatCellPair_2(i);
    end
end
    
%% plot all cells after remove the repeat cells
%
close all
CMP=jet(5);
figure
width=1400;
height=1400;
x0=10;
y0=10;
set(gcf,'position',[x0,y0,width,height])
a3=gca;
k=1;
Plane=2;
CMP=[0,1,0;1,0,1];
for h=1:1:size(Cell_P2,1)
        FOV_ID=NeuronMatrix(Cell_P2(h),1);
        i=NeuronMatrix(Cell_P2(h),2);
        if ~ismember(h,RepeatCell)
                Xpoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.xpix)',256);
                Ypoint=mod(double(StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellStat{1,i}.ypix)',256);
                Xpoint(Xpoint==0)=256;
                Ypoint(Ypoint==0)=256;
                ROI_pannel=zeros(size(Image_mean));
                for m=1:1:length(Xpoint)
                    ROI_pannel(Ypoint(m),Xpoint(m))=1;
                end
                ROI_pannel_corrected=imwarp(ROI_pannel,TransformMatrix,'OutputView',imref2d(size(ROI_pannel)));             
                ROI_pannel_corrected= imrotate(ROI_pannel_corrected,90);
                ROI_pannel_corrected = flip(ROI_pannel_corrected,2);
                if FOV_ID==2
                ROI_pannel_corrected = imrotate(ROI_pannel_corrected,0.6); 
                else
                end
                [Ypoint,Xpoint,~]=find(ROI_pannel_corrected>0);      
                P=convhull(Xpoint,Ypoint);             
                fill(FOVposition_correct(FOV_ID,1)+Xpoint(P),FOVposition_correct(FOV_ID,2)+Ypoint(P),CMP(FOV_ID,:),'LineStyle','none','facealpha',.6);                    
                hold on;    
            else
        end
       h
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
%     camroll(-90)
    set(gca,'YDir','reverse');
    set(gcf,'color',[0 0 0]);
    set(gca,'color',[0 0 0]);
    set(a3, 'XAxisLocation', 'top')
    set(a3, 'TickDir', 'out')
box off
axis off
