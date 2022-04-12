classdef MINI2P < dabs.resources.Device & most.HasMachineDataFile & dabs.resources.configuration.HasConfigPage & dabs.resources.widget.HasWidget
    properties (SetAccess=protected,Hidden)
        ConfigPageClass = 'dabs.resources.configuration.resourcePages.MINI2PPage';
    end
    
    methods (Static)
        function names = getDescriptiveNames()
            names = {'MINI2P'};
        end
    end
    
    properties (SetAccess=protected)
        WidgetClass = 'dabs.resources.widget.widgets.MINI2PWidget'; 
    end
    
    %% ABSTRACT PROPERTY REALIZATIONS (most.HasMachineDataFile)
    properties (Constant, Hidden)
        %Value-Required properties
        mdfClassName = mfilename('class');
        mdfHeading = 'MINI2P Settings';
        
        %Value-Optional properties
        mdfDependsOnClasses; %#ok<MCCPI>
        mdfDirectProp;       %#ok<MCCPI>
        mdfPropPrefix;       %#ok<MCCPI>
        
        mdfDefault = defaultMdfSection();
    end
    
    properties (SetObservable)
        transformMatrix;
        transformMatrixDirectory;
        scope;
        system;
        objective;
    end
    
    methods
        function obj = MINI2P(name)
            obj@dabs.resources.Device(name);
            obj@most.HasMachineDataFile(true);
            
            obj.deinit();
            obj.loadMdf();
            obj.reinit();
        end
        
        function delete(obj)
            obj.deinit();
        end
    end
    
    methods
        function reinit(obj)
            obj.deinit();
            try
                %Figure out how to install matlab app if not already
                %installed
                A=matlab.apputil.getInstalledAppInfo;
                B=0;
                ID=[];
                for i=1:1:size(A,2)
                    if strcmp('MINI2PDistortionDetection',A(1,i).name)
                        B=1;
                        ID=i;
                    else
                    end
                end
                
                if B==0
                    [path,~,~]=fileparts(which('scanimage'))
                    matlab.apputil.install...
                        ([path,'\plugin\MINI2PDistortionDetection.mlappinstall']);
                    disp('MINI2PDistortionDetection is installed now!');
                else
%                     matlab.apputil.uninstall(A(1,i).id);
%                     disp('Previous MINI2PDistortionDetection is removed!');
%                     matlab.apputil.install...
%                         ('C:\Program Files\Vidrio\SI2021.0.0_2021-10-06-114206_3baa4b74be\+dabs\+mini2p\MINI2PDistortionDetection.mlappinstall');
                    disp('MINI2PDistortionDetection is already installed!');
                    
                end
    
                obj.errorMsg = ''; 
            catch ME
                obj.deinit();
                obj.errorMsg = sprintf('%s: initialization error: %s',obj.name,ME.message);
                most.ErrorHandler.logError(ME,obj.errorMsg);
            end   
        end
        
        function deinit(obj)            
           % Put device back into a cleared state

            obj.errorMsg = 'uninitialized';
        end
    end
    
    methods        
      
        function loadMdf(obj)
            success = true;
            success = success & obj.safeSetPropFromMdf('transformMatrixDirectory', 'transformMatrixDirectory');
            success = success & obj.safeSetPropFromMdf('system', 'system');            
            success = success & obj.safeSetPropFromMdf('scope', 'scope');
            success = success & obj.safeSetPropFromMdf('objective', 'objective');
            if ~success
                obj.errorMsg = 'Error loading config';
            end
        end
        
        function saveMdf(obj)
            
            obj.safeWriteVarToHeading('transformMatrixDirectory', obj.transformMatrixDirectory);
            obj.safeWriteVarToHeading('system', obj.system);
            obj.safeWriteVarToHeading('scope', obj.scope);
            obj.safeWriteVarToHeading('objective', obj.objective);
        end
    end
    
    %% User Methods - Core functionality of the device
    methods
        function distortionCorrection (obj, varargin)
            % select data to process
            hSI = evalin('base', 'hSI');
            hSICtl = evalin('base', 'hSICtl');
            Recentpath=hSI.hScan2D.logFilePath;
            DataFolder = uigetdir(Recentpath);
            Tifffiles=dir([DataFolder,'\*.tif']);
            for i=1:1:size(Tifffiles,1)
                disp(Tifffiles(i,1).name);
            end
            ZoomD=str2double(hSICtl.hGUIs.mainControlsV4.Children(24,1).Children(9,1).String)/10;
            ZoomOnes=str2double(hSICtl.hGUIs.mainControlsV4.Children(24,1).Children(11,1).String);
            ZoomTens=str2double(hSICtl.hGUIs.mainControlsV4.Children(24,1).Children(12,1).String)*10;
            Zoom=ZoomD+ZoomOnes+ZoomTens;
            % register 3D tranformation matrix path
            if isempty(obj.transformMatrixDirectory)
                [file,path] = uigetfile('*.*','Select Transform Matrix File');
                
                if isnumeric(file)
                    return % user cancelled
                end
                file = fullfile(path,file);
                TranformMatrix=load(file);
                
                obj.transformMatrixDirectory = file;
                obj.saveMdf();
            else
                TranformMatrix = load(obj.transformMatrixDirectory);
            end
            % start to do distortion correction
            f = waitbar(0,'Start to load data...');
            h=1;
            tStart = tic;
            % Aout: MxNxCxFxSxV array, with images of size MxN for C channels, F frames, S slices, and V volumes. Default type is uint16.
            for p=1:1:size(Tifffiles,1)
                %     waitbar(h/imgInfo.numSlices,f,['Tranforming:',num2str(h),'/',num2str(imgInfo.numSlices),' (',num2str(T,'%.0f'),' s)']);
                disp(['Start load: "' Tifffiles(p,1).name]);
                [header, RawImage, imgInfo] = scanimage.util.opentif([Tifffiles(p,1).folder,'\',Tifffiles(p,1).name]);
                disp(['MINI2P imaging data: "' Tifffiles(p,1).name  '" is loaded']);
                % first check if this is corrected data or raw data;
                ms=[];
                if isfield(header.SI, 'Custom')
                    if header.SI.Custom.MINI2P.MINI2P_Corrected==1
                        disp([Tifffiles(p,1).name,' is already corrected. Jump to the next file!']);
                        ms=msgbox([Tifffiles(p,1).name,' is already corrected.Jump to the next file!']);
%                         drawnow;
%                         pause(10)
                       
                    else
                        disp([Tifffiles(p,1).name,' was applied an unknown modification. Correction ignored!']);
                        ms=msgbox([Tifffiles(p,1).name,' was applied an unknown modification. Correction ignored!']);
%                         drawnow;
%                         pause(10)
                    end
                else
                    ImagingSize=size(RawImage);
                    Zs=header.SI.hStackManager.zsRelative;
                    Zs=unique(Zs,'stable');
                    % check which transformax to use (nearest Z)
                    Zdistance=zeros(size(TranformMatrix.TransfMatrix,1),imgInfo.numSlices);
                    for i=1:1:size(TranformMatrix.TransfMatrix,1)
                        Zdistance(i,:)=Zs-TranformMatrix.TransfMatrix{i,1}.Z;
                    end                   
                    ID=[];
                    for i=1:1:imgInfo.numSlices
                        [~,ID(i)]=min(abs(Zdistance(:,i)));
                    end
                    % correct plan one by one
                    
                    
                    TransformMatrix_used=zeros(imgInfo.numSlices,1);                                       
                    for i=1:1:imgInfo.numSlices
                        TransformMatrix_tem=TranformMatrix.TransfMatrix{ID(i),1};
                        TransformMatrix_used(i,1)=ID(i);
                        if Zoom~=TransformMatrix_tem.Zoom
                            answer = questdlg('Current Zoom is not equal to the Zoom saved in transform Maxtix!', ...
                                'Warning', ...
                                'Continue','Cancel','Cancel');
                            if strcmp (answer, 'Cancel')
                                disp('Correction is abort');
                                return
                            else
                            end
                        else
                        end
                        % check if use the original matrax, or need to scale it to fit the image size
                        Scale=size(RawImage,1)/str2double(TransformMatrix_tem.LineNum);
                        if Scale~=1
                            Anchor_pairs_withdistortion=TransformMatrix_tem.GridParameter.Anchor_pairs_withdistortion*Scale;
                            Anchor_pairs_nodistortion=TransformMatrix_tem.GridParameter.Anchor_pairs_nodistortion*Scale;
                            FOV=TransformMatrix_tem.FOVx;
                            Res=TransformMatrix_tem.Resx/Scale;
                            TransformMatrix_final=fitgeotrans(Anchor_pairs_withdistortion,Anchor_pairs_nodistortion,'pwl');
                        else
                            TransformMatrix_final=TransformMatrix_tem.Trasfm;
                        end
                        T=toc(tStart);
                        waitbar(h/((imgInfo.numSlices*size(Tifffiles,1))+1),f,['Transforming: Tiff ',num2str(p),'/ Plane ',num2str(i),' (',num2str(T,'%.0f'),' s)']);
                        RawImage(:,:,:,:,i,:)=imwarp( RawImage(:,:,:,:,i,:),TransformMatrix_final(1),'OutputView',imref2d(size( RawImage(:,:,:,:,i,:))));
                        h=h+1;
                    end
                    
                    
                    cd (Tifffiles(p,1).folder)
                    disp(['Saving...',Tifffiles(p,1).name(1:end-4),'_wrappiing corrected.tif']);
                    
                    
                    % also put the transform matrix and the correction tag
                    % into tif header.
                    
                    customHdr=[sprintf('SI.Custom.MINI2P.MINI2P_Corrected = true\n') ...
                        sprintf('SI.Custom.MINI2P.tf_path = ''%s''\n',obj.transformMatrixDirectory) ...
                        sprintf('SI.Custom.MINI2P.tfused = ''%s''\n',num2str(TransformMatrix_used')) ...
                        sprintf('SI.Custom.MINI2P.MINI2P_system = ''%s''\n',TransformMatrix_tem.system) ...
                        sprintf('SI.Custom.MINI2P.MINI2P_scope = ''%s''\n',TransformMatrix_tem.scope) ...
                        sprintf('SI.Custom.MINI2P.MINI2P_objective = ''%s''\n',TransformMatrix_tem.objective)];
                    
                    
                                       
                    resaveTiff2(Tifffiles(p,1).name, RawImage, [Tifffiles(p,1).name(1:end-4),'_distcorrected.tif'],customHdr)
                    
                    
                end
%                 delete (ms);
            end
            T=toc(tStart);
            waitbar(1,f,['Transform finished!']);
            
            
            function resaveTiff2(fileName, AoutImdata, newFileName, customHdr)
                if nargin < 4 || isempty(customHdr)
                    customHdr = [];
                end
                
                hTif = scanimage.util.ScanImageTiffReader(fileName);
                [fileHeaderStr,frameDescs] = getHeaderDataFromScanImageTiffObj(hTif);
                hTif.close();
                
                [header,~,~,~] = scanimage.util.opentif(fileName);
                [hMroiRoiGroup,~,~] = scanimage.util.readTiffRoiData(fileName);
                
                %     v = extractBetween(fileHeaderStr, 'framesPerSlice', sprintf('\n'));
                %     val = ['SI.hStackManager.framesPerSlice' v{1} sprintf('\n')];
                %
                %     [startIndex,endIndex] = regexp(fileHeaderStr, val);
                %
                %     b4 = fileHeaderStr(1:startIndex-1);
                %     after = fileHeaderStr(endIndex:end);
                %
                %     new = [b4 sprintf('SI.hStackManager.framesPerSlice = %d', size(AoutImdata,4)) after];
                
                %%
                [startIdx, ~] = regexp(fileHeaderStr, 'SI.hStackManager.framesPerSlice');
                b4 = fileHeaderStr(1:startIdx-1);
                endIdx = regexp(fileHeaderStr, 'SI.hStackManager.name');
                after = fileHeaderStr(endIdx:end);
                
                new = [b4 sprintf('SI.hStackManager.framesPerSlice = %d\n', size(AoutImdata,4)) after];
                
                fileHeaderStr = new;
                
                [startIdx, ~] = regexp(fileHeaderStr, 'SI.hStackManager.numSlices');
                b4 = fileHeaderStr(1:startIdx-1);
                endIdx = regexp(fileHeaderStr, 'SI.hStackManager.numVolumes');
                after = fileHeaderStr(endIdx:end);
                
                new = [b4 sprintf('SI.hStackManager.numSlices = %d\n', size(AoutImdata,5)) after];
                
                fileHeaderStr = new;
                
                [startIdx, ~] = regexp(fileHeaderStr, 'SI.hStackManager.numVolumes');
                b4 = fileHeaderStr(1:startIdx-1);
                endIdx = regexp(fileHeaderStr, 'SI.hStackManager.reserverInfo');
                after = fileHeaderStr(endIdx:end);
                
                new = [b4 sprintf('SI.hStackManager.numVolumes = %d\n', size(AoutImdata,6)) after];
                
                fileHeaderStr = [customHdr new];
                
                %%
                
                hdrSplit = strsplit(fileHeaderStr, '\n\n');
                
                hdrBuf = [uint8(hdrSplit{1}) 0];
                hdrBufLen = length(hdrBuf);
                
                hdrBuf = [hdrBuf uint8([sprintf('\n') hdrSplit{2}]) 0];
                
                pfix = [1 3 3 7 typecast(uint32(4),'uint8') typecast(uint32(hdrBufLen),'uint8') typecast(uint32(length(hdrBuf)-hdrBufLen),'uint8')];
                
                tifHeaderData = [pfix hdrBuf]';
                tifHeaderStringOffset = length(pfix);
                tifRoiDataStringOffset = length(pfix) + hdrBufLen;
                
                
                dataSigned    = true;
                bitsPerSample = 16;
                
                numChannelSave = numel(header.SI.hChannels.channelSave);
                pixelsPerLine = header.SI.hRoiManager.pixelsPerLine;
                linesPerFrame = header.SI.hRoiManager.linesPerFrame;
                
                blankFrameDescription = repmat(' ',1,2000);
                
                imageSize = pixelsPerLine * linesPerFrame * (bitsPerSample/8);
                
                flybackPeriods = ceil(header.SI.hScan2D.flybackTimePerFrame * header.SI.hScan2D.scannerFrequency);
                flybackLinesPerFrame = round((flybackPeriods * 2^header.SI.hScan2D.bidirectional)/2)*2;
                
                rois = hMroiRoiGroup.rois;
                
                zs = header.SI.hStackManager.zs;
                
                scanFields = arrayfun(@(z)hMroiRoiGroup.scanFieldsAtZ(z),...
                    zs,'UniformOutput',false);
                
                cumPixelResolutionAtZ = zeros(0,2);
                mRoiLogging = false;
                for zidx = 1:length(scanFields)
                    sfs = scanFields{zidx};
                    pxRes = zeros(0,2);
                    for sfidx = 1:length(sfs)
                        sf = sfs{sfidx};
                        pxRes(end+1,:) = sf.pixelResolution(:)';
                    end
                    mRoiLogging = mRoiLogging || size(pxRes,1) > 1;
                    cumPixelResolutionAtZ(end+1,:) = [max(pxRes(:,1)),sum(pxRes(:,2))+((size(pxRes,1)-1)*flybackLinesPerFrame)];
                end
                
                mRoiLogging = mRoiLogging || any(cumPixelResolutionAtZ(1,1) ~= cumPixelResolutionAtZ(:,1));
                mRoiLogging = mRoiLogging || any(cumPixelResolutionAtZ(1,2) ~= cumPixelResolutionAtZ(:,2));
                linesPerFrame = max(cumPixelResolutionAtZ(:,2));
                pixelsPerLine = max(cumPixelResolutionAtZ(:,1));
                
                sf = scanFields{1}{1};
                resDenoms = 2^30 ./ (1e4 * sf.pixelResolutionXY ./ (sf.sizeXY * header.SI.objectiveResolution));
                
                xResolutionNumerator = 2^30;
                xResolutionDenominator = resDenoms(1);
                yResolutionNumerator = 2^30;
                yResolutionDenominator = resDenoms(2);
                
                imageSize = pixelsPerLine * linesPerFrame * (bitsPerSample/8);
                
                hNewTif = scanimage.components.scan2d.TiffStream;
                
                assert(hNewTif.open([pwd '\' newFileName],tifHeaderData,tifHeaderStringOffset,tifRoiDataStringOffset), 'Failed to create log file.');
                hNewTif.configureImage(pixelsPerLine, linesPerFrame, (bitsPerSample/8), numChannelSave, dataSigned, blankFrameDescription,...
                    xResolutionNumerator, xResolutionDenominator, yResolutionNumerator, yResolutionDenominator);
                
                frameStream = {};
                
                %‘Aout’ is specified, a matrix of the size M x N x C x F x S x V is
                % created, where: - C spans the channel indices, - F spans the frame
                % indicies, - S spans the slice indices, and - V the volume indices.
                
                numChans = size(AoutImdata,3);
                numFrames = size(AoutImdata,4);
                numSlices = size(AoutImdata,5);
                numVolumes = size(AoutImdata,6);
                
                for vols = 1:numVolumes
                    for slice = 1:numSlices
                        for frames = 1:numFrames
                            for chan = 1:numChans
                                frameStream{end+1} = AoutImdata(:,:,chan, frames, slice, vols)';
                            end
                        end
                    end
                end
                
                
                for idx = 1:numel(frameStream)
                    
                    hNewTif.replaceImageDescription(frameDescs{idx});
                    
                    if mRoiLogging
                        line = 1;
                        for roiIdx = 1:length(rois)
                            imdata = frameStream{idx};
                            dims = size(imdata);
                            tempbuf(1:dims(1),line:line+dims(2)-1) = imdata;
                            line = line + dims(2);
                        end
                        hNewTif.appendFrame(int16(tempbuf), imageSize);
                    else
                        hNewTif.appendFrame(int16(frameStream{idx}), imageSize);
                    end
                end
                
                hNewTif.close();
                hNewTif.cleanUp();
                
                function [fileHeaderStr,frameDescs] = getHeaderDataFromScanImageTiffObj(tifObj)
                    frameDescs = tifObj.descriptions();
                    isemptyMask = cellfun(@(d)isempty(d),frameDescs);
                    frameDescs(isemptyMask) = [];
                    
                    fileHeaderStr = tifObj.metadata();
                end
            end


            
            
            
        end
        
    end
    
    %% Setters and Getters Validation
    methods
        
    end
end



function s = defaultMdfSection()
s = [...
    most.HasMachineDataFile.makeEntry('transformMatrixDirectory' ,'','transform Matrix Directory terminal  e.g. ''C:\Users\''')...
    most.HasMachineDataFile.makeEntry('system' ,'','Name of the system  e.g. ''System_001''')...
    most.HasMachineDataFile.makeEntry('scope' ,'','Name of the scope  e.g. ''MINI2P-L_001''')...
    most.HasMachineDataFile.makeEntry('objective' ,'','Type of the objective  e.g. ''D0233''')...
    ];
end




% ----------------------------------------------------------------------------
% Version1: (2022.03.22) at Kavli Institute, Trondheim, NTNU
% The MINI2P distortion detection and correction method is from the paper Zong, et.al,Cell,2022
% The MINI2P distortion detection GUI was written by Weijian Zong in the Moser lab
% The SI device of MINI2P was modified from standard device in SI by Weijian Zong in the Moser lab, Mitchell Sandoe and Jacob Franklin in Vidrio Technologies, LLC
% Updated version of this device and the user manual can be found in https://github.com/kavli-ntnu/MINI2P_toolbox
% ----------------------------------------------------------------------------
% Copyright (C) 2021 Vidrio Technologies, LLC
% 
% ScanImage (R) 2021 is software to be used under the purchased terms
% Code may be modified, but not redistributed without the permission
% of Vidrio Technologies, LLC
% 
% VIDRIO TECHNOLOGIES, LLC MAKES NO WARRANTIES, EXPRESS OR IMPLIED, WITH
% RESPECT TO THIS PRODUCT, AND EXPRESSLY DISCLAIMS ANY WARRANTY OF
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
% IN NO CASE SHALL VIDRIO TECHNOLOGIES, LLC BE LIABLE TO ANYONE FOR ANY
% CONSEQUENTIAL OR INCIDENTAL DAMAGES, EXPRESS OR IMPLIED, OR UPON ANY OTHER
% BASIS OF LIABILITY WHATSOEVER, EVEN IF THE LOSS OR DAMAGE IS CAUSED BY
% VIDRIO TECHNOLOGIES, LLC'S OWN NEGLIGENCE OR FAULT.
% CONSEQUENTLY, VIDRIO TECHNOLOGIES, LLC SHALL HAVE NO LIABILITY FOR ANY
% PERSONAL INJURY, PROPERTY DAMAGE OR OTHER LOSS BASED ON THE USE OF THE
% PRODUCT IN COMBINATION WITH OR INTEGRATED INTO ANY OTHER INSTRUMENT OR
% DEVICE.  HOWEVER, IF VIDRIO TECHNOLOGIES, LLC IS HELD LIABLE, WHETHER
% DIRECTLY OR INDIRECTLY, FOR ANY LOSS OR DAMAGE ARISING, REGARDLESS OF CAUSE
% OR ORIGIN, VIDRIO TECHNOLOGIES, LLC's MAXIMUM LIABILITY SHALL NOT IN ANY
% CASE EXCEED THE PURCHASE PRICE OF THE PRODUCT WHICH SHALL BE THE COMPLETE
% AND EXCLUSIVE REMEDY AGAINST VIDRIO TECHNOLOGIES, LLC.
% ----------------------------------------------------------------------------
