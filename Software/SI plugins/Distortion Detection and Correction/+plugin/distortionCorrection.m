function distortionCorrection ()
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
TranformMatrix=load('C:\Recordings\Distortion correction test\file_00001_TransformMatrix_MINI2P-L_001_1.0_0307223552.mat');
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
    for i=1:1:imgInfo.numSlices
        TransformMatrix_tem=TranformMatrix.TransfMatrix{ID(i),1};
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
    plugin.resaveTiff2(Tifffiles(p,1).name, RawImage, [Tifffiles(p,1).name(1:end-4),'_wrappiing corrected.tif'])
end    
T=toc(tStart);
waitbar(1,f,['Transform finished!']); 

function resaveTiff2(fileName, AoutImdata, newFileName)
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
    
    fileHeaderStr = new;
    
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
