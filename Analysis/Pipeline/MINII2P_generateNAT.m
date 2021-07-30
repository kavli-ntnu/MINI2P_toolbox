%% this is the first part of the analysis of MINI2P imaging data.
close all
clear all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2P data processing%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% step 1: set all paramters for reading and processing the data 
RaveSIInformation=0; % if read and save the scanimage information to the same forder.
SaveRaw=0; % if save raw imaging data into the output dataset. Choose SaveRaw=0, if you want to save space and processing time.
SaveNeuronActiveMatrix=1;% if save NeuronActiveMatrix. 
SaveTrackingResultMatrix=1;% if save TrackingResultMatrix
SaveNAT=1; % if save FullNeuronBehaviorDataSet.
ReadTrackingInforfromMetaFile=1; % if directly read tracking information (box size, box center,etc) fromt the metafile.
Input_CalciumLowpassFreqency=2; % Unit:Hz, the cutting frenquency for the lowpass filter for calium signal;
Input_MovingWindow=30;  % Unit:s, the time window for calculating the local baseline signal (F-zero). The origin of this number was taken from paper:"Cellular resolution optical access to brain regions in fissures Imaging medial prefrontal cortex and grid cells in entorhinal cortex" PNAS,(2014)
Input_TrancientTimethreadhold=0.75; % Unit:s, how long do the significant trancients need to last. The purpose for setting this parameter is to filter out the fake trancients by motion artifacts and the detector noise.  
Input_timeexpanded=0.4; % Unit:s, give a small time before and after the detected significant trancient. The purpose of this time is to catch the complete trancient.
Input_CellDistance_threadhold=10; % unit: pixel, distance threadhold for clastering the same cell in two planes;
Input_CellCrossCorrelation_threadhold=0.5; % unit: a.u., cross-correlation threadhold for clastering the same cell in two planes;
Input_Overlap_threadhold=0.8;% unit: % the overlap region/cell RIO
Input_TrackingBodyPart={};
Input_TrackingBodyPart{1,1}={'nose','mouse','lefthand','righthand','leftleg','rightleg','tailbase','bodycenter'}; % the body parts tracked by DLC bottom camera
Input_TrackingBodyPart{2,1}={'miniscope','leftear','rightear','bodycenter','tailbase'}; % the body parts tracked by DLC top camera
Input_TrackingBodyPart{3,1}={'nose','leftear','rightear','bodycenter','tailbase'}; % the body parts tracked by DLC top camera in Emilie setup
Input_NeuronpilFactor=0.7;
Input_multisessioncombined=1; % if the data is going to be analyzed as a combined session
Input_TrackingFrameratesource=1; %1 for sourse from avi video; 2 for sourse from tracking metafile, 3 for from the suite2p output, 4 for munually input;
Input_TrackingCorrection=4;

% bellow input value are only used if no metafile is provided 
Input_ImagingFramesize_default=[256 256]; %2P imaging frame size
Input_FocalPlane_default=[0 -40];
Input_cameratype_default=2;
Input_BoxScale_default=0.9; % cm/pixel
Input_BoxSize_default=80; % cm
Input_BoxCenter_default=[558,614]; % middle of the frame, 800x800
Input_framerate_manually=15; % tracking frame rate
%%
%%%%%%%%%%% step2: select the raw file root folder %%%%%%%%%%%

File_folder_raw=uigetdir; % munually select the raw file folder;

if File_folder_raw~=0
    
    disp('Data folder was founded');
else
    disp('Data folder does not exist');
end
cd (File_folder_raw); % the working folder will be changed to the data folder
%
%%%%%%%%%%%% step3: load in suite2P output %%%%%%%%%%%%

% if Multiplane imaging was applied, the combined folder should be
% generated. Then only the data in combined folder will be analyzed.
% Otherwise, only analyze the plane0 data.
if exist([File_folder_raw,'\suite2p\combined'],'dir')==7  
    Suite2PFolder=[File_folder_raw,'\suite2p\combined\'];
    [NeuronInformation,NeuronActiveMatrix]=preprocessing.GenerateNeuronInformation([File_folder_raw,'\suite2p\combined\Fall.mat']);
else
    Suite2PFolder=[File_folder_raw,'\suite2p\plane0\'];
    [NeuronInformation,NeuronActiveMatrix]=preprocessing.GenerateNeuronInformation([File_folder_raw,'\suite2p\plane0\Fall.mat']);
end
% save the folder information to the NeuronInformation package
NeuronInformation.File_folder_Suite2P=Suite2PFolder;
NeuronInformation.File_folder_root=File_folder_raw;
%%
%clear up the non-cell based on Suite2P classifier
% if TrustSuite2PIdentifier==1
%     for i=1:1:NeuronInformation.TotalROINumber
%         if NeuronInformation.CellPossibility(i,2)==0
  
%%%%%%%%%% step4 (optional): load in the unwrapping transform Matric %%%%%
%%%%%%%%%% 
% load('C:\PROGRAMS\MINICalium_beta\MINI2P_analysis\ReferenceDATA\UnwrappingMatrix\OPENMINI2000Hz-001-Zoom1-256x256.mat')
% NeuronInformation.TransformMatrix=TransformMatrix;
% step5 (optional): save the SI information to seperate files
%header has the detail imaging information including the hardware setups and recoding conditions
%Aout{i} is the raw 16bit data, organized by the size M x N x C x F x S x V , where C spans the channel indices, F spans the frame indicies, S spans the slice indices, and V the volume indices. 
%imgInfor includes the general imaging infortmation
% search for all .tif files in this folder, only necessasy for SI information saving. No need for suite2p or DLC output analysis
NeuronInformation.TIF_Files=dir([NeuronInformation.File_folder_root,'\*.tif']);
[TIF_num ~]=size(NeuronInformation.TIF_Files);
if TIF_num~=0
    disp([num2str(TIF_num),' Imaging tif files were found']);
    for i=1:1:TIF_num
        disp(NeuronInformation.TIF_Files(i).name);
    end
else
    disp('0 imaging tif files were not found');
end
NeuronActiveMatrix.RawImage={}; %the content for the raw imaging data
NeuronInformation.header=[];
NeuronInformation.imgInfo=[];

if RaveSIInformation==1
    if SaveRaw==1
        for WhichTiff=1:1:size(NeuronInformation.TIF_Files,1)
            disp(['Now loading file "' num2str(WhichTiff)  '.....']);
            [header, NeuronInformation.RawImage{WhichTiff}, imgInfo] = scanimage.util.opentif([NeuronInformation.TIF_Files(WhichTiff).folder,'\',NeuronInformation.TIF_Files(WhichTiff).name]);
            % [header1, Aout1, imgInfo1] = scanimage.util.opentif([TIF_Files(WhichTiff).folder,'\',TIF_Files(1).name]);
            filename1=imgInfo.filename;
            NeuronInformation.header=header;
            NeuronInformation.imgInfo=imgInfo;
            disp(['2P imaging data: "' imgInfo  '" is loaded']);
            save([filename1(1:end-4),'_tifHeader.mat'],'header');% save the SI information, only necessasy for SI information saving. No need for DJ analysis
            disp([filename1(1:end-4),'_tifHeader.mat',' is saved']);
            save([filename1(1:end-4),'_imgInfo.mat'],'imgInfo');% save the imaging information,only necessasy for SI information saving. No need for DJ analysis
            disp([filename1(1:end-4),'_imgInfo.mat',' is saved']);
        end
    else
            for WhichTiff=1:1:size(NeuronInformation.TIF_Files,1)
            disp(['Now loading file "' num2str(WhichTiff)  '.....']);
            [header,~, imgInfo] = scanimage.util.opentif([NeuronInformation.TIF_Files(WhichTiff).folder,'\',NeuronInformation.TIF_Files(WhichTiff).name]);
            NeuronInformation.header=header;
            NeuronInformation.imgInfo=imgInfo;            
            filename1=imgInfo.filename;
            disp(['2P imaging data: "' filename1  '" is loaded']);
            save([filename1(1:end-4),'_tifHeader.mat'],'header');% save the SI information, only necessasy for SI information saving. No need for DJ analysis
            disp([filename1(1:end-4),'_tifHeader.mat',' is saved']);
            save([filename1(1:end-4),'_imgInfo.mat'],'imgInfo');% save the imaging information,only necessasy for SI information saving. No need for DJ analysis
            disp([filename1(1:end-4),'_imgInfo.mat',' is saved']);
            end
    end
else
end

clear header
clear imgInfo1
clear filename1
clear WhichTiff
clear TIF_num

%%%%%%%%%%%% step6: pre-processing the calcium data %%%%%%%%%%
% first tick out the cells idenfied as "Not-cell" in Suite2p
NeuronInformation.TotalCell=length(find(NeuronInformation.CellPossibility(:,1)==1));
NeuronInformation.IsCell=find(NeuronInformation.CellPossibility(:,1)==1);
NeuronInformation.ROIStat_Iscell=NeuronInformation.ROIStat(1,NeuronInformation.IsCell);
NeuronInformation.ROIStat_Iscell=NeuronInformation.ROIStat(1,NeuronInformation.IsCell);
NeuronActiveMatrix.EventTrain_IsCell=zeros(size(NeuronActiveMatrix.EventTrain,1),NeuronInformation.TotalCell+2);
NeuronActiveMatrix.EventTrain_IsCell(:,1:2)=NeuronActiveMatrix.EventTrain(:,1:2);
NeuronActiveMatrix.EventTrain_IsCell(:,3:NeuronInformation.TotalCell+2)=NeuronActiveMatrix.EventTrain(:,NeuronInformation.IsCell+2);
NeuronActiveMatrix.F_raw_IsCell=zeros(size(NeuronActiveMatrix.F_raw,1),NeuronInformation.TotalCell+2);
NeuronActiveMatrix.F_raw_IsCell(:,1:2)=NeuronActiveMatrix.F_raw(:,1:2);
NeuronActiveMatrix.F_raw_IsCell(:,3:NeuronInformation.TotalCell+2)=NeuronActiveMatrix.F_raw(:,NeuronInformation.IsCell+2);
NeuronActiveMatrix.F_neuropil_IsCell=zeros(size(NeuronActiveMatrix.F_neuropil,1),NeuronInformation.TotalCell+2);
NeuronActiveMatrix.F_neuropil_IsCell(:,1:2)=NeuronActiveMatrix.F_neuropil(:,1:2);
NeuronActiveMatrix.F_neuropil_IsCell(:,3:NeuronInformation.TotalCell+2)=NeuronActiveMatrix.F_neuropil(:,NeuronInformation.IsCell+2);

% calculate the raw detaF/F, threadholded deconvolved spike 
    NeuronInformation.LowpassFreqency=Input_CalciumLowpassFreqency; %Hz
    NeuronInformation.MovingWindow=Input_MovingWindow; %second
    NeuronInformation.TrancientTimethreadhold=Input_TrancientTimethreadhold; %second
    NeuronInformation.timeexpanded=Input_timeexpanded;
    NeuronInformation.NeuronpilFactor=Input_NeuronpilFactor;
    NeuronActiveMatrix.detaF=NeuronActiveMatrix.F_raw_IsCell;
    NeuronActiveMatrix.detaF_F=NeuronActiveMatrix.F_raw_IsCell;
    NeuronActiveMatrix.detaF_F_filtered=NeuronActiveMatrix.F_raw_IsCell;
    NeuronActiveMatrix.F_zero=NeuronActiveMatrix.F_raw_IsCell;
    NeuronActiveMatrix.SignificantTrancient=NeuronActiveMatrix.F_raw_IsCell;
    NeuronActiveMatrix.detaF(:,3:end)=NeuronActiveMatrix.F_raw_IsCell(:,3:end)-Input_NeuronpilFactor*NeuronActiveMatrix.F_neuropil_IsCell(:,3:end); % neuropil correction   
%     NeuronActiveMatrix.detaF(NeuronActiveMatrix.detaF(:)<0)=0;
    if ~isempty(NeuronActiveMatrix.RawImage)
       NeuronInformation.Framesize=size(NeuronActiveMatrix.RawImage{1},[1 2]); % suppose all images were taken with the size frame size
    else
       NeuronInformation.Framesize=Input_ImagingFramesize_default;  % if raw data was not saved in the NeuronActiveMatrix, then the frame size will be manually input
    end

%calculate the focal plane information for each cell 
    NeuronInformation.CellinPlane=zeros(NeuronInformation.TotalCell,1); 
    NeuronInformation.FocalPlane=zeros(NeuronInformation.ImagingPlane,1); 
    if ~isempty(NeuronInformation.header) % get the focal plane from the SI meta data         
        if NeuronInformation.ImagingPlane>1
            NeuronInformation.FocalPlane=NeuronInformation.header.SI.hStackManager.zsRelative;
        else
            NeuronInformation.FocalPlane=0; % right now just manually set the focal plane if only one plane was recorded. Should be able to read this directly from SI metadata;
        end  
    else    % set the focal plane manually
        NeuronInformation.FocalPlane=Input_FocalPlane_default;
    end

%Calcium trace pre-processing 
    for i=3:1:size(NeuronActiveMatrix.detaF_F,2)
        [NeuronActiveMatrix.detaF_F_filtered(:,i) NeuronActiveMatrix.detaF_F(:,i) NeuronActiveMatrix.F_zero(:,i) NeuronActiveMatrix.SignificantTrancient(:,i)]=preprocessing.DetaF_T_baselinecorrect(NeuronActiveMatrix.detaF(:,i),NeuronInformation.LowpassFreqency,NeuronInformation.MovingWindow,NeuronInformation.TrancientTimethreadhold,NeuronInformation.volumerate,NeuronInformation.timeexpanded);
        disp(['Baseline correction for neuron: ', num2str(i)]);
    end    
    NeuronActiveMatrix.detaF_F_filtered=NeuronActiveMatrix.detaF_F;   
    for i=3:1:NeuronInformation.TotalCell+2
        NeuronActiveMatrix.detaF_F_filtered(:,i)=NeuronActiveMatrix.detaF_F(:,i).*NeuronActiveMatrix.SignificantTrancient(:,i); % filtered calcium trace by picking up the significnat trancient.
    end  
% deconvolved calcium event  filteting and scaling to df/f 
    NeuronActiveMatrix.FilteredEvent(:,1:2)=NeuronActiveMatrix.detaF_F(:,1:2);    
    
    for i=3:1:NeuronInformation.TotalCell+2
        Event_baseline=find(NeuronActiveMatrix.SignificantTrancient(:,i)==0);     
        Event_uncertainty=mean(NeuronActiveMatrix.EventTrain_IsCell(Event_baseline,i));
        NeuronActiveMatrix.FilteredEvent(:,i)=(NeuronActiveMatrix.EventTrain_IsCell(:,i)-Event_uncertainty).*NeuronActiveMatrix.SignificantTrancient(:,i); % only keep the deconvoled spikes within the significant trancient period.
        for j=1:1:size(NeuronActiveMatrix.FilteredEvent(:,i),1)
            if NeuronActiveMatrix.FilteredEvent(j,i)<=0
                NeuronActiveMatrix.FilteredEvent(j,i)=0;
            else
            end
        end  
        % remove outlier
        Outliers=find(NeuronActiveMatrix.FilteredEvent(:,i)>2*prctile(NeuronActiveMatrix.FilteredEvent(:,i),99.9));    
        NeuronActiveMatrix.FilteredEvent(Outliers,i)=prctile(NeuronActiveMatrix.FilteredEvent(:,i),99.95);    
        dF_F_MAX=max(NeuronActiveMatrix.detaF_F_filtered(:,i));

        Event_MAX=max(NeuronActiveMatrix.FilteredEvent(:,i));
        NeuronActiveMatrix.FilteredEvent(:,i)=NeuronActiveMatrix.FilteredEvent(:,i).*dF_F_MAX./Event_MAX;  % normalized the Event to df/f, important step!!!!!
    end
    
% clearup bad frames (No signal)
NeuronInformation.BadCell=[];
% 
k=1;
for i=1:1:NeuronInformation.TotalCell
    if sum(NeuronActiveMatrix.F_raw_IsCell(:,i+2))==0
        NeuronInformation.BadCell(k)=i;
        k=k+1;
    else
    end
end

% Calculate the SNR of each cell
    % the SNR is defined as the average of peak amptitutes (the first peak after each deconvoled Event) devides the average STD in non-trancient period of the calcium trace; 
    NeuronInformation.SNR=zeros(NeuronInformation.TotalCell,1);    
    for i=1:1:NeuronInformation.TotalCell
        NeuronInformation.SNR(i)=preprocessing.CalculateSNR(NeuronActiveMatrix.detaF_F(:,i+2),find(NeuronActiveMatrix.SignificantTrancient(:,i+2)==1),find(NeuronActiveMatrix.SignificantTrancient(:,i+2)==0));
%         NeuronInformation.SNR(i)=preprocessing.CalculateSNR(NeuronActiveMatrix.detaF_F(:,i+2),find(NeuronActiveMatrix.SignificantTrancient(:,i+2)==1),find(NeuronActiveMatrix.SignificantTrancient(:,i+2)==0));        
        disp(['SNR of neuron ',num2str(i),' is ',num2str(NeuronInformation.SNR(i))]);
        pause(0.05);
    end
% calculate the Event Count
NeuronInformation.EventCount_raw=zeros(NeuronInformation.TotalCell,1);
%
    for i=1:1:NeuronInformation.TotalCell
        NeuronInformation.EventCount_raw(i,1)=length(find(NeuronActiveMatrix.FilteredEvent(:,i+2)>0));
        disp(['Events count of neuron ',num2str(i),' is ',num2str(NeuronInformation.EventCount_raw(i,1))]);
    end
%
%%%%%% claster the same cell in different planes %%%%%%

% the method of deciding the same cell in different plane is：
% Step1：calculat the cross-correlation of each cell pairs in two adjecent
% plane；
% Step2：calculate the center distance of each cell pairs in two adjecent
% plane；
% Step3：calculate the overlap ratio of the cell footprints
% Step4: set the threadholds of cross-correlation, center distance and overlapping;
% Step5: link the threadholded cell pairs to the same cell family. The
% famility tree can cross multiplanes;
% Step6: keep only the brightest cells in a family; 

for i=1:1:NeuronInformation.TotalCell
    NeuronInformation.ROIStat_Iscell{1,i}.med(1)=mod(NeuronInformation.ROIStat_Iscell{1,i}.med(1),NeuronInformation.Framesize(1)); % recover the real cell position from the combined (tailored) data (data in combined folder)
    NeuronInformation.ROIStat_Iscell{1,i}.med(2)=mod(NeuronInformation.ROIStat_Iscell{1,i}.med(2),NeuronInformation.Framesize(2)); % recover the real cell position from the combined (tailored) data (data in combined folder)
end

if NeuronInformation.ImagingPlane>1
    for i=1:1:NeuronInformation.TotalCell
        NeuronInformation.CellinPlane(i,1)=NeuronInformation.ROIStat_Iscell{i}.iplane;
    end
    NeuronInformation.CellinEachPlane=zeros(NeuronInformation.ImagingPlane,1);
    for i=1:1:NeuronInformation.ImagingPlane
        NeuronInformation.CellinEachPlane(i)=size(find(NeuronInformation.CellinPlane(:,1)==(i-1)),1);
    end
else
    for i=1:1:NeuronInformation.TotalROINumber
        NeuronInformation.CellinPlane(i,1)=0;
    end
    NeuronInformation.CellinEachPlane=zeros(NeuronInformation.ImagingPlane,1);
    for i=1:1:NeuronInformation.ImagingPlane
        NeuronInformation.CellinEachPlane(i)=size(find(NeuronInformation.CellinPlane(:,1)==(i-1)),1);
    end
end

NeuronInformation.CellCrosscorrelation_betweenPlane=cell(1,NeuronInformation.ImagingPlane-1);
NeuronInformation.CellCenter_betweenPlane=cell(1,NeuronInformation.ImagingPlane-1);
NeuronInformation.CellOverlap_betweenPlane=cell(1,NeuronInformation.ImagingPlane-1);
NeuronInformation.StartCell_eachPlane=preprocessing.CalculateStartFrame(NeuronInformation.CellinEachPlane', 1);
if NeuronInformation.ImagingPlane>1
    for i=1:1:NeuronInformation.ImagingPlane-1
        NeuronInformation.CellCrosscorrelation_betweenPlane{1,i}=zeros(NeuronInformation.CellinEachPlane(i),NeuronInformation.CellinEachPlane(i+1));
        NeuronInformation.CellCenter_betweenPlane{1,i}=zeros(NeuronInformation.CellinEachPlane(i),NeuronInformation.CellinEachPlane(i+1));        
        for j=1:1:NeuronInformation.CellinEachPlane(i)
            for k=1:1:NeuronInformation.CellinEachPlane(i+1)
                r=corr(smoothdata(NeuronActiveMatrix.F_raw_IsCell(:,NeuronInformation.StartCell_eachPlane(i)+j+1),'gaussian',3),smoothdata(NeuronActiveMatrix.F_raw_IsCell(:,NeuronInformation.StartCell_eachPlane(i+1)+k+1),'gaussian',3));
                NeuronInformation.CellCrosscorrelation_betweenPlane{1,i}(j,k)=r;
                NeuronInformation.CellCenter_betweenPlane{1,i}(j,k)=sqrt(sum((NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i)+j-1}.med-NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i+1)+k-1}.med).^2));
                OverLay=preprocessing.overlap([NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i)+j-1}.ypix;NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i)+j-1}.xpix],[NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i+1)+k-1}.ypix;NeuronInformation.ROIStat_Iscell{1,NeuronInformation.StartCell_eachPlane(i+1)+k-1}.xpix],NeuronInformation.Framesize);
                NeuronInformation.CellOverlap_betweenPlane{1,i}(j,k)=max(OverLay(1:2));
           end
            disp(['search for child cell for cell number ',num2str(j),' in plane ',num2str(i)]);
        end
    end
else
end
%
NeuronInformation.Overlap_threadhold=Input_Overlap_threadhold;
NeuronInformation.CellDistance_threadhold=Input_CellDistance_threadhold;
NeuronInformation.CellCrossCorrelation_threadhold=Input_CellCrossCorrelation_threadhold;
SameCellCandidate=cell(1,NeuronInformation.ImagingPlane-1);

for i=1:1:NeuronInformation.ImagingPlane-1   
    [x y]=find((NeuronInformation.CellCenter_betweenPlane{1,i}(:,:)<NeuronInformation.CellDistance_threadhold).*(NeuronInformation.CellCrosscorrelation_betweenPlane{1,i}(:,:)>NeuronInformation.CellCrossCorrelation_threadhold).*(NeuronInformation.CellOverlap_betweenPlane{1,i}(:,:)>NeuronInformation.Overlap_threadhold));
    NeuronInformation.SameCellCandidate{1,i}=[x y];
end

NeuronInformation.Samecelltracking=zeros(1,NeuronInformation.TotalCell);
% 
for i=1:1:NeuronInformation.ImagingPlane-1
    for j=1:1:size(NeuronInformation.SameCellCandidate{1,i},1)
        NeuronInformation.Samecelltracking(NeuronInformation.StartCell_eachPlane(i)+NeuronInformation.SameCellCandidate{1,i}(j,1)-1)=NeuronInformation.StartCell_eachPlane(i+1)+NeuronInformation.SameCellCandidate{1,i}(j,2)-1;
    end
end

NeuronInformation.SamecellFamily={};
family_number=1;
Samecelltracking_temp=NeuronInformation.Samecelltracking;
% familytree=find(Samecelltracking_temp~=0);
% SamecellFamily_tem=zeros(2,size(familytree,2));
% link cells to the same famaly tree in more than 2 planes
for i=1:1:size(Samecelltracking_temp,2)
    if Samecelltracking_temp(i)~=0
        Mum=i;
        child=Samecelltracking_temp(Mum);
        family=[Mum,child];
        while Samecelltracking_temp(child)~=0
        family=[family,Samecelltracking_temp(child)];    
        child=Samecelltracking_temp(child);
        end
        NeuronInformation.SamecellFamily{1,family_number}=family;
        family_number=family_number+1;
        Samecelltracking_temp(family)=0;
    else
    end
end
%
NeuronInformation.SamecellFamily_rawID=cell(size(NeuronInformation.SamecellFamily));
NeuronInformation.SamecellFamily_Brightest=zeros(size(NeuronInformation.SamecellFamily));
NeuronInformation.RepeatCell=[];
for i=1:1:size(NeuronInformation.SamecellFamily,2)
    for j=1:1:size(NeuronInformation.SamecellFamily{1,i},2)
        NeuronInformation.SamecellFamily_rawID{1,i}(j)=NeuronInformation.IsCell(NeuronInformation.SamecellFamily{1,i}(j));
        BrightnessinFamily(j)=mean(NeuronActiveMatrix.detaF_F(find(NeuronActiveMatrix.SignificantTrancient(:,NeuronInformation.SamecellFamily{1,i}(j)+2)>0),NeuronInformation.SamecellFamily{1,i}(j)+2));   
    end
    [~,BrightestID]=max(BrightnessinFamily);
    NeuronInformation.SamecellFamily_Brightest(i)=NeuronInformation.SamecellFamily{1,i}(BrightestID);
    NeuronInformation.RepeatCell=[NeuronInformation.RepeatCell,setdiff(NeuronInformation.SamecellFamily{1,i},NeuronInformation.SamecellFamily{1,i}(BrightestID))];
    
end
NeuronInformation.RepeatCell=unique(NeuronInformation.RepeatCell);    
% show calcium activity, SNR and events counts for all cells
close all
CellPerPlot=60;
Offset=-5;% the inteval between each traces, unit: df/f
k=1;
Xrange=[0 500];
% for j=1:CellPerPlot:NeuronInformation.TotalCell
for j=1:CellPerPlot:NeuronInformation.TotalCell
    figure (k)    
    x0=10+(k-1)*420;
    y0=100;
    width=400;
    height=1200;
    set(gcf,'position',[x0,y0,width,height]) 
    m=1;
    for i=j:1:min(j+CellPerPlot-1,NeuronInformation.TotalCell)
        ID=i;
        plot(NeuronActiveMatrix.F_raw_IsCell(:,1),smooth((m-1)*Offset+NeuronActiveMatrix.detaF_F(:,ID+2),5),'color',[0 0 0],'Linewidth',1);
        hold on
        Significant=find(NeuronActiveMatrix.SignificantTrancient(:,ID+2)==1);
        scatter(NeuronActiveMatrix.detaF_F(Significant,1),(m-1)*Offset+NeuronActiveMatrix.detaF_F(Significant,ID+2),2*ones(size(Significant)),[1 0 0],'filled')
        text(Xrange(1)-80,(m-1)*Offset,[num2str(ID),' (',num2str(NeuronInformation.IsCell(ID)),')'],'FontSize',7);
        text(Xrange(2)+10,(m-1)*Offset,num2str(NeuronInformation.SNR(ID),'%.1f'),'FontSize',7);
        m=m+1;
    end
ylim([CellPerPlot*Offset 5])
xlim([Xrange])
box off
axis off
k=k+1;
title(['Cell ',num2str(j),' to Cell ',num2str(min(j+CellPerPlot-1,NeuronInformation.TotalCell))]);
end

    figure (k)
    x0=10+(k-1)*420;
    y0=100;
    width=1200;
    height=800;
    subplot(2,1,1)
    set(gcf,'position',[x0,y0,width,height]) 
    plot(NeuronInformation.SNR)
    title('SNR for each cell');
    subplot(2,1,2)
    plot(NeuronInformation.EventCount_raw)
    title('Event count for each cell');
%%
%%%%%% output NeuronInformation package and NeuronActivy matrix  %%%%%%
%
if SaveNeuronActiveMatrix==1
    save ([NeuronInformation.File_folder_root,'\NeuronActiveMatrix.mat'],'NeuronActiveMatrix');
    save ([NeuronInformation.File_folder_root,'\NeuronInformation.mat'],'NeuronInformation');
else
    disp('NueronInformation and NeuronAcitivity files are not saved!!!!!');
end
disp('pre-processing of the neuronal activity data was done!');
%%
%%%%%%%%%%%%%%%%%%%% process the tarcking data %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step1 set the tracking parameters from the raw tracking video, and load in the tracking output
TrackingResultMatrix=struct;
TrackingResultMatrix.TrackingVideoPath=dir([File_folder_raw,'\*trackingVideo.avi']);
TrackingResultMatrix.SessionRAW=size(TrackingResultMatrix.TrackingVideoPath,1);
TrackingResultMatrix.File_folder_raw=File_folder_raw;
for i=1:1:TrackingResultMatrix.SessionRAW
    Video=VideoReader([TrackingResultMatrix.TrackingVideoPath(i,1).folder,'/',TrackingResultMatrix.TrackingVideoPath(i,1).name]);
    TrackingResultMatrix.FrameineachSessionRAW(i)=Video.NumFrames;
    TrackingResultMatrix.VideoWidthineachSessionRAW(i)=Video.Width;
    TrackingResultMatrix.VideoHeightineachSessionRAW(i)=Video.Height;
end
%
if Input_multisessioncombined==1  
    TrackingResultMatrix.Session=1;% combined all session    
    TrackingResultMatrix.FrameineachSession=sum(TrackingResultMatrix.FrameineachSessionRAW);
    disp('all sessions are combined to one session for analysis!');
else
    TrackingResultMatrix.Session=TrackingResultMatrix.SessionRAW; %extract from suip2p output
    TrackingResultMatrix.FrameineachSession=TrackingResultMatrix.FrameineachSessionRAW;%extract from suip2p output
    disp('each session are analyzing seperately!');
end
%
%%%%%% get the tracking enviroment information %%%%%%
TrackingResultMatrix.TrackingMetaFile=dir([File_folder_raw,'\*MetaFile.csv']);
if ReadTrackingInforfromMetaFile==1
    Meta=tdfread(TrackingResultMatrix.TrackingMetaFile(1).name);
    TrackingResultMatrix.cameratype=str2num(Meta.Value(1,:));  
    TrackingResultMatrix.BoxScale=str2num(Meta.Value(7,:));% cm/pixel     
    TrackingResultMatrix.BoxSize=str2num(Meta.Value(8,:)); % cm
    TrackingResultMatrix.BoxCenter(1)=str2num(Meta.Value(9,:)); %pixel, middle of the box
    TrackingResultMatrix.BoxCenter(2)=str2num(Meta.Value(10,:)); %pixel, middle of the box
%     TrackingResultMatrix.cameratype=str2num(Meta.Value0x22(1,1:find(Meta.Value0x22(1,:)=='"')-1));  
%     TrackingResultMatrix.BoxScale=str2num(Meta.Value0x22(7,1:find(Meta.Value0x22(7,:)=='"')-1));% cm/pixel     
%     TrackingResultMatrix.BoxSize=str2num(Meta.Value0x22(8,1:find(Meta.Value0x22(8,:)=='"')-1)); % cm
%     TrackingResultMatrix.BoxCenter(1)=str2num(Meta.Value0x22(9,1:find(Meta.Value0x22(9,:)=='"')-1)); %pixel, middle of the box
%     TrackingResultMatrix.BoxCenter(2)=str2num(Meta.Value0x22(10,1:find(Meta.Value0x22(10,:)=='"')-1)); %pixel, middle of the box
else
    TrackingResultMatrix.cameratype=Input_cameratype_default;
    TrackingResultMatrix.BoxScale=Input_BoxScale_default; % cm/pixel
    TrackingResultMatrix.BoxSize=Input_BoxSize_default; % cm
    TrackingResultMatrix.BoxCenter=Input_BoxCenter_default; % middle of the frame, 800x800
end
%
if Input_TrackingFrameratesource==1
    TrackingResultMatrix.framerate=Video.FrameRate;%Hz,get the frame rate from avi video
elseif Input_TrackingFrameratesource==2        
    TrackingResultMatrix.framerate=str2num(Meta.Value0x22(4,[1:find(BBB.Value0x22(4,:)=='"')-1])); % Hz, get the frame rate from tracking metafile
elseif Input_TrackingFrameratesource==3
    TrackingResultMatrix.framrate=double(NeuronActiveMatrix.Options.fs)*double(ImagingPlane);  %Hz, get the frame rate from the suite2p output
elseif Input_TrackingFrameratesource==4
    TrackingResultMatrix.framrate=Input_framerate_manually; % mannully input
else
end
%
TrackingResultMatrix.BodypartRaw={};
TrackingResultMatrix.BodypartRaw=Input_TrackingBodyPart{TrackingResultMatrix.cameratype,1};
TrackingResultMatrix.TimeStamp=cell(1,TrackingResultMatrix.Session);
TrackingResultMatrix.DLCfilename=cell(1,TrackingResultMatrix.Session);

for i=1:1:TrackingResultMatrix.Session
    TrackingResultMatrix.TimeStamp{1,i}=double((1:1:TrackingResultMatrix.FrameineachSession(i))-1)./TrackingResultMatrix.framerate;
    TrackingResultMatrix.TimeStamp{1,i}=TrackingResultMatrix.TimeStamp{1,i}';
end
%
%%%%%% load in the DLC output %%%%%%s
TrackingResultMatrix.DLCRaw=cell(1,TrackingResultMatrix.SessionRAW);
for i=1:1:TrackingResultMatrix.SessionRAW
    disp(['please select DLC output file ', num2str(i)]);
    TrackingResultMatrix.DLCfilename{1,i}=uigetfile('*.csv');
    disp(['load in DLC output: ', TrackingResultMatrix.DLCfilename{1,i}]);
    TrackingResultMatrix.DLCRaw{1,i}=csvread(TrackingResultMatrix.DLCfilename{1,i},3,0);
end

DLCCombined=[];
TrackingResultMatrix.DLC={};
if Input_multisessioncombined==1
    for i=1:1:TrackingResultMatrix.SessionRAW
        DLCCombined=[DLCCombined;TrackingResultMatrix.DLCRaw{1,i}];
    end
    TrackingResultMatrix.DLC{1,1}=DLCCombined;
else
    TrackingResultMatrix.DLC=TrackingResultMatrix.DLCRaw;
end
%
%%%%%% clear up the raw data by appling a filter for highly cenfitatial rate %%%%%%
TrackingResultMatrix.DLCFiltered=cell(1,TrackingResultMatrix.Session);
TrackingResultMatrix.DLCFiltered_percent=cell(1,TrackingResultMatrix.Session);
for i=1:1:TrackingResultMatrix.Session
   [ TrackingResultMatrix.DLCFiltered{1,i},TrackingResultMatrix.DLCFiltered_percent{1,i}]=preprocessing.adp_filt(TrackingResultMatrix.DLC{1,i});
    
end
%
TrackingResultMatrix.TrackingdataRaw=cell(1,TrackingResultMatrix.Session); %%seperate diffirent body parts from cleaned up DLC data,box size callibrated, and centered
if TrackingResultMatrix.cameratype==1
    for i=1:1:TrackingResultMatrix.Session
        TrackingdataRaw_tem=TrackingResultMatrix.DLCFiltered{1,i};
        TrackingdataRaw_tem(:,1:2:end)=(TrackingdataRaw_tem(:,1:2:end)-TrackingResultMatrix.BoxCenter(1))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingdataRaw_tem(:,2:2:end)=-1*(TrackingdataRaw_tem(:,2:2:end)-TrackingResultMatrix.BoxCenter(2))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingResultMatrix.TrackingdataRaw{1,i}.noseposition=TrackingdataRaw_tem(:,1:2);
        TrackingResultMatrix.TrackingdataRaw{1,i}.mouseposition=TrackingdataRaw_tem(:,3:4);
        TrackingResultMatrix.TrackingdataRaw{1,i}.lefthandposition=TrackingdataRaw_tem(:,5:6);
        TrackingResultMatrix.TrackingdataRaw{1,i}.righthandposition=TrackingdataRaw_tem(:,7:8);
        TrackingResultMatrix.TrackingdataRaw{1,i}.leftlegposition=TrackingdataRaw_tem(:,9:10);   
        TrackingResultMatrix.TrackingdataRaw{1,i}.rightlegposition=TrackingdataRaw_tem(:,11:12);    
        TrackingResultMatrix.TrackingdataRaw{1,i}.tailbaseposition=TrackingdataRaw_tem(:,13:14);   
        TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition=TrackingdataRaw_tem(:,15:16);    
    end
elseif TrackingResultMatrix.cameratype==2
    for i=1:1:TrackingResultMatrix.Session
        TrackingdataRaw_tem=TrackingResultMatrix.DLCFiltered{1,i};
        TrackingdataRaw_tem(:,1:2:end)=(TrackingdataRaw_tem(:,1:2:end)-TrackingResultMatrix.BoxCenter(1))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingdataRaw_tem(:,2:2:end)=-1*(TrackingdataRaw_tem(:,2:2:end)-TrackingResultMatrix.BoxCenter(2))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingResultMatrix.TrackingdataRaw{1,i}.miniscopeposition=TrackingdataRaw_tem(:,1:2);
        TrackingResultMatrix.TrackingdataRaw{1,i}.leftearposition=TrackingdataRaw_tem(:,3:4);
        TrackingResultMatrix.TrackingdataRaw{1,i}.rightearposition=TrackingdataRaw_tem(:,5:6);
        TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition=TrackingdataRaw_tem(:,7:8);
        TrackingResultMatrix.TrackingdataRaw{1,i}.tailbaseposition=TrackingdataRaw_tem(:,9:10);     
    end
    
elseif TrackingResultMatrix.cameratype==3
    for i=1:1:TrackingResultMatrix.Session
        TrackingdataRaw_tem=TrackingResultMatrix.DLCFiltered{1,i};
        TrackingdataRaw_tem(:,1:2:end)=(TrackingdataRaw_tem(:,1:2:end)-TrackingResultMatrix.BoxCenter(1))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingdataRaw_tem(:,2:2:end)=-1*(TrackingdataRaw_tem(:,2:2:end)-TrackingResultMatrix.BoxCenter(2))*TrackingResultMatrix.BoxScale; %change unit from pixel to cm, and center the box
        TrackingResultMatrix.TrackingdataRaw{1,i}.noseposition=TrackingdataRaw_tem(:,15:16);
        TrackingResultMatrix.TrackingdataRaw{1,i}.leftearposition=TrackingdataRaw_tem(:,9:10);
        TrackingResultMatrix.TrackingdataRaw{1,i}.rightearposition=TrackingdataRaw_tem(:,11:12);
        TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition=(TrackingdataRaw_tem(:,15:16)+TrackingdataRaw_tem(:,13:14))./2;
        TrackingResultMatrix.TrackingdataRaw{1,i}.tailbaseposition=TrackingdataRaw_tem(:,13:14);     
    end    
    
  
end
%

%%%%%% get the head and body position and direction %%%%%%
TrackingResultMatrix.PositionDirectionExtract=cell(1,TrackingResultMatrix.Session);

for i=1:1:TrackingResultMatrix.Session
TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition=zeros(TrackingResultMatrix.FrameineachSession(i),2);
TrackingResultMatrix.PositionDirectionExtract{1,i}.headdirection=zeros(TrackingResultMatrix.FrameineachSession(i),1);
TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyposition=zeros(TrackingResultMatrix.FrameineachSession(i),2);
TrackingResultMatrix.PositionDirectionExtract{1,i}.bodydirection=zeros(TrackingResultMatrix.FrameineachSession(i),1);
TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodyangle=zeros(TrackingResultMatrix.FrameineachSession(i),1);
TrackingResultMatrix.PositionDirectionExtract{1,i}.headspeed=zeros(TrackingResultMatrix.FrameineachSession(i),1);
TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyspeed=zeros(TrackingResultMatrix.FrameineachSession(i),1);
end
%
if TrackingResultMatrix.cameratype==1
% calculate the headposition
% the center of nose and mouse,centered in box middle
    for i=1:1:TrackingResultMatrix.Session
        TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition=(TrackingResultMatrix.TrackingdataRaw{1,i}.noseposition(:,1:2)+TrackingResultMatrix.TrackingdataRaw{1,i}.mouseposition(:,1:2))/2;
    end
% calculate the headdirection
% the line from mouse to nose, east 0,north 90, 
    for i=1:1:TrackingResultMatrix.Session
        for j=1:1:TrackingResultMatrix.FrameineachSession(i)
            TrackingResultMatrix.PositionDirectionExtract{1,i}.headdirection(j)=SpatialTuning_BNT.DirectionCalulate(TrackingResultMatrix.TrackingdataRaw{1,i}.mouseposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.mouseposition(j,2),TrackingResultMatrix.TrackingdataRaw{1,i}.noseposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.noseposition(j,2));
        end
    end
else
 % calculate the headposition
% the center of leftear and rightear,centered in box middle
    for i=1:1:TrackingResultMatrix.Session
        TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition=(TrackingResultMatrix.TrackingdataRaw{1,i}.leftearposition(:,1:2)+TrackingResultMatrix.TrackingdataRaw{1,i}.rightearposition(:,1:2))/2;
    end
% calculate the headdirection
% the perpendicular line to the line from leftear to right( mouse to nose), east 0,north 90, 
    for i=1:1:TrackingResultMatrix.Session
        for j=1:1:TrackingResultMatrix.FrameineachSession(i)
            TrackingResultMatrix.PositionDirectionExtract{1,i}.headdirection(j)=mod(SpatialTuning_BNT.DirectionCalulate(TrackingResultMatrix.TrackingdataRaw{1,i}.leftearposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.leftearposition(j,2),TrackingResultMatrix.TrackingdataRaw{1,i}.rightearposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.rightearposition(j,2))+90,360);
        end
    end   
    
end
% calculate the bodyposition
% the position of bodycenter;
for i=1:1:TrackingResultMatrix.Session
    TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyposition=TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition(:,1:2);
end
% calculate the bodydirection
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session
    for j=1:1:TrackingResultMatrix.FrameineachSession(i)
    TrackingResultMatrix.PositionDirectionExtract{1,i}.bodydirection(j)=SpatialTuning_BNT.DirectionCalulate(TrackingResultMatrix.TrackingdataRaw{1,i}.tailbaseposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.tailbaseposition(j,2),TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition(j,1),TrackingResultMatrix.TrackingdataRaw{1,i}.bodycenterposition(j,2));
    end
end

% calculate the headbodyangle
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session   
    TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodyangle=TrackingResultMatrix.PositionDirectionExtract{1,i}.bodydirection-TrackingResultMatrix.PositionDirectionExtract{1,i}.headdirection;
    TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodyangle=rem(TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodyangle,180);
end
% calculate the headspeed
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session  
    TrackingResultMatrix.PositionDirectionExtract{1,i}.headspeed=SpatialTuning_BNT.SpeedCalculate(TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition,1./TrackingResultMatrix.framerate,15);
end
% calculate the bodyspeed
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session  
    TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyspeed=SpatialTuning_BNT.SpeedCalculate(TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyposition,1./TrackingResultMatrix.framerate,15);
end
% calculate the headbodyspeed
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session  
    TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodyspeed=TrackingResultMatrix.PositionDirectionExtract{1,i}.headspeed-TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyspeed;
end
% calculate the headbodydistance
% the line from tailbase to bodycenter, east 0,north 90, 
for i=1:1:TrackingResultMatrix.Session
    TrackingResultMatrix.PositionDirectionExtract{1,i}.headbodydistance=BNTanalysis_modified.DistanceCalculate(TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition,TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyposition);
end
%

%%%%%% output dataset %%%%%%
%
if SaveTrackingResultMatrix==1
        save ([TrackingResultMatrix.File_folder_raw,'\TrackingResultMatrix.mat'],'TrackingResultMatrix');     
else
    disp('tracking data analysis was not saved!!!!')
end
disp('tracking data analysis was done!');

%% generate NAT

NAT=struct;
ExperimentInformation=struct;
AAA=find(NeuronInformation.File_folder_root=='\');
ExperimentInformation.Animal=NeuronInformation.File_folder_root(AAA(end-2)+1:AAA(end-1)-1);
ExperimentInformation.Experiment=NeuronInformation.File_folder_root(AAA(end-1)+1:AAA(end)-1);
ExperimentInformation.Date=NeuronInformation.File_folder_root(AAA(end)+1:end);
ExperimentInformation.Session=TrackingResultMatrix.Session; % how many session in one recording
ExperimentInformation.BoxScale=TrackingResultMatrix.BoxScale;  % mm/pixel
ExperimentInformation.BoxSize=TrackingResultMatrix.BoxSize;  % mm/pixel
ExperimentInformation.BoxCenter=TrackingResultMatrix.BoxCenter; 
ExperimentInformation.RawDataAddress=NeuronInformation.File_folder_root;
ExperimentInformation.FrameRate=TrackingResultMatrix.framerate;
ExperimentInformation.FrameinEachSession=TrackingResultMatrix.FrameineachSession;
ExperimentInformation.TimeStamp=TrackingResultMatrix.TimeStamp;
ExperimentInformation.DLCfilename=TrackingResultMatrix.DLCfilename;
ExperimentInformation.cameratype=TrackingResultMatrix.cameratype;
ExperimentInformation.BodypartRaw=TrackingResultMatrix.BodypartRaw;
ExperimentInformation.ImagingChannel=NeuronInformation.ImagingChannel;
ExperimentInformation.ImagingPlane=NeuronInformation.ImagingPlane;
ExperimentInformation.TotalCell=NeuronInformation.TotalCell;
ExperimentInformation.IsCell=NeuronInformation.IsCell;
ExperimentInformation.CellinPlane=NeuronInformation.CellinPlane;
ExperimentInformation.CellinEachPlan=NeuronInformation.CellinEachPlane;
ExperimentInformation.CellSNR=NeuronInformation.SNR;
ExperimentInformation.EventCount_raw=NeuronInformation.EventCount_raw;
ExperimentInformation.CellStat=NeuronInformation.ROIStat_Iscell;
ExperimentInformation.FocalPlane=NeuronInformation.FocalPlane;
ExperimentInformation.Framesize=NeuronInformation.Framesize;
ExperimentInformation.FrameShift=Input_TrackingCorrection;
ExperimentInformation.ImagingOptions=NeuronInformation.Options;
ExperimentInformation.RepeatCell=NeuronInformation.RepeatCell;
NAT=cell(1,ExperimentInformation.Session);


%
%colume 1:  timestamp
%colume 2:3 headposition
%colume 4:  headdirection
%colume 5:  headspeed
%colume 6:  headvilid
%colume 7:8 bodyposition
%colume 9:  bodydirection
%colume 10:  bodyspeed
%colume 11:  bodyvilid
%colume 12: behavoir clastering
%colume 13:4:13+N*4-1 neuron detaF/F_raw, detaF/F_filtered, detaF_F_significant and deconvoled spikes
ConfidenceThreadhold=0.1;
Startframe=preprocessing.CalculateStartFrame(ExperimentInformation.FrameinEachSession,double(ExperimentInformation.ImagingPlane));
for i=1:1:ExperimentInformation.Session
    NAT{1,i}=NaN(ExperimentInformation.FrameinEachSession(i),13+ExperimentInformation.TotalCell*4-1);
    NAT{1,i}(:,1)=ExperimentInformation.TimeStamp{i};
    NAT{1,i}(:,2:3)=TrackingResultMatrix.PositionDirectionExtract{1,i}.headposition;

    NAT{1,i}(:,4)=TrackingResultMatrix.PositionDirectionExtract{1,i}.headdirection;
    NAT{1,i}(:,5)=TrackingResultMatrix.PositionDirectionExtract{1,i}.headspeed;   
    if ExperimentInformation.cameratype==1 
        for j=1:1:ExperimentInformation.FrameinEachSession(i)
            if (TrackingResultMatrix.DLC{1,i}(j,4)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,7)<ConfidenceThreadhold)
                NAT{1,i}(j,6)=0;
            else
                NAT{1,i}(j,6)=1;
            end
            if (TrackingResultMatrix.DLC{1,i}(j,22)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,25)<ConfidenceThreadhold)
                NAT{1,i}(j,11)=0;
            else
                NAT{1,i}(j,11)=1;
            end
        end
    elseif ExperimentInformation.cameratype==2 
       for j=1:1:ExperimentInformation.FrameinEachSession(i)
            if (TrackingResultMatrix.DLC{1,i}(j,7)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,10)<ConfidenceThreadhold)
                NAT{1,i}(j,6)=0;
            else
                NAT{1,i}(j,6)=1;
            end
            if (TrackingResultMatrix.DLC{1,i}(j,13)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,15)<ConfidenceThreadhold)
                NAT{1,i}(j,11)=0;
            else
                NAT{1,i}(j,11)=1;
            end
        end
    elseif ExperimentInformation.cameratype==3
        for j=1:1:ExperimentInformation.FrameinEachSession(i)
                    if (TrackingResultMatrix.DLC{1,i}(j,16)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,19)<ConfidenceThreadhold)
                NAT{1,i}(j,6)=0;
            else
                NAT{1,i}(j,6)=1;
            end
            if (TrackingResultMatrix.DLC{1,i}(j,22)<ConfidenceThreadhold)||(TrackingResultMatrix.DLC{1,i}(j,25)<ConfidenceThreadhold)
                NAT{1,i}(j,11)=0;
            else
                NAT{1,i}(j,11)=1;
            end
        end
    else
    end
    NAT{1,i}(:,2)=circshift(NAT{1,i}(:,2),-1*ExperimentInformation.FrameShift);
    NAT{1,i}(:,3)=circshift(NAT{1,i}(:,3),-1*ExperimentInformation.FrameShift);  
    NAT{1,i}(:,7:8)=TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyposition;
    NAT{1,i}(:,9)=TrackingResultMatrix.PositionDirectionExtract{1,i}.bodydirection;
    NAT{1,i}(:,10)=TrackingResultMatrix.PositionDirectionExtract{1,i}.bodyspeed;
    NAT{1,i}(:,7)=circshift(NAT{1,i}(:,7),-1*ExperimentInformation.FrameShift);
    NAT{1,i}(:,8)=circshift(NAT{1,i}(:,8),-1*ExperimentInformation.FrameShift);
        if i==ExperimentInformation.Session
            for j=1:1:ExperimentInformation.TotalCell 
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,9+j*4)=NeuronActiveMatrix.detaF_F(Startframe(i):end,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,10+j*4)=NeuronActiveMatrix.detaF_F_filtered(Startframe(i):end,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,11+j*4)=NeuronActiveMatrix.SignificantTrancient(Startframe(i):end,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,12+j*4)=NeuronActiveMatrix.FilteredEvent(Startframe(i):end,j+2);
            end
        else
            for j=1:1:ExperimentInformation.TotalCell             
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,9+j*4)=NeuronActiveMatrix.detaF_F(Startframe(i):Startframe(i+1)-1,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,10+j*4)=NeuronActiveMatrix.detaF_F_filtered(Startframe(i):Startframe(i+1)-1,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,11+j*4)=NeuronActiveMatrix.SignificantTrancient(Startframe(i):Startframe(i+1)-1,j+2);
                NAT{1,i}(ExperimentInformation.CellinPlane(j)+1:ExperimentInformation.ImagingPlane:end,12+j*4)=NeuronActiveMatrix.FilteredEvent(Startframe(i):Startframe(i+1)-1,j+2);
            end
        end
        
end 
%

%%%%%% output dataset %%%%%%
if SaveNAT==1
        save ([ExperimentInformation.RawDataAddress,'\NAT.mat'],'NAT','-v7.3');
        save ([ExperimentInformation.RawDataAddress,'\ExperimentInformation.mat'],'ExperimentInformation','-v7.3');
else
end
disp('NAT generation was done!');

close all



