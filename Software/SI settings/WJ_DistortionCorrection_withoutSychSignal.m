function WJ_DistortionCorrection_withoutSychSignal(src,evt,varargin)
            global filename        
            filename=[src.hSI.hScan2D.logFilePath,'\', src.hSI.hScan2D.logFileStem,'_',num2str(src.hSI.hScan2D.logFileCounter-1,'%05d'),'.tif'];                                   
            disp(['Just saved raw data: ',filename]);
            disp('Start the distortion correction....');
            disp('Reloading the raw data....');
            % read Tranform matrix
            load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\MINI2P\FOV callibration\50umGrid_256_TransformMatrix.mat');
            % read FOV information
            FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\MINI2P\FOV callibration\50umGrid-256x256-2000Hz_FOV.csv');
            % Open the .tif recording file desired to apply the distortion correction
            [header, RawImage, imgInfo] = scanimage.util.opentif('C:\Recordings\Experiment1\256_data.tif');
            filename1=imgInfo.filename;
            disp(['Rawe 2P imaging data: "' filename1  '" is loaded']);
            save([filename(1:end-4),'_tifHeader.mat'],'header');% save the SI information, only necessasy for SI information saving. No need for DJ analysis
            disp([filename(1:end-4),'_tifHeader.mat',' is saved']);
            save([filename(1:end-4),'_imgInfo.mat'],'imgInfo');% save the imaging information,only necessasy for SI information saving. No need for DJ analysis
            disp([filename(1:end-4),'_imgInfo.mat',' is saved']);
            Dimention=size(RawImage); 
            if length(Dimention)<=3                  
                ImageStack_corrected=RawImage;
                disp('Image Stack Size:');
                disp(num2str(Dimention));
                disp('No HyperStack Merging is applied');
            elseif length(Dimention)==4 
                disp('Image Stack Size:');
                disp(num2str(Dimention));
                disp('2D Timelapse imaging was infered.  Channel and Frame merging is applied');
                ImageStack_corrected=zeros(Dimention(1),Dimention(2),Dimention(3)*Dimention(4));
                k=1;
                for j=1:1:Dimention(4)
                    for i=1:1:Dimention(3)
                        ImageStack_corrected(:,:,k)=RawImage(:,:,i,j);   
                        k=k+1;
                    end
                end
            elseif length(Dimention)==5
                disp('Image Stack Size:');
                disp(num2str(Dimention));
                disp('Multi-layer sinlge TimePoint imaging was infered.  Channel, Frame and Plane merging is applied');
                ImageStack_corrected=zeros(Dimention(1),Dimention(2),Dimention(3)*Dimention(4)*Dimention(5));    
                k=1;
                    for m=1:1:Dimention(5)
                        for j=1:1:Dimention(4)
                            for i=1:1:Dimention(3)
                                ImageStack_corrected(:,:,k)=RawImage(:,:,i,j,m);   
                                k=k+1;
                            end
                        end
                    end
            else
                disp('Image Stack Size:');
                disp(num2str(Dimention));
                disp('Multi-layer Timelapse imaging was infered.  Channel, Frame, Plane and Volume merging is applied');
                ImageStack_corrected=zeros(Dimention(1),Dimention(2),Dimention(3)*Dimention(4)*Dimention(5)*Dimention(6));
                k=1;
                for n=1:1:Dimention(6)
                    for m=1:1:Dimention(5)
                        for j=1:1:Dimention(4)
                            for i=1:1:Dimention(3)
                                ImageStack_corrected(:,:,k)=RawImage(:,:,i,j,m,n);   
                                k=k+1;
                            end
                        end
                    end
                end    
            end
            clear RawImage
            fTIF = DataIO.Fast_BigTiff_Write([filename(1:end-4),'_wrappiing-corrected.tif'],1*10000/FOV(3),0);
            for i=1:1:size(ImageStack_corrected,3)               
                ImageStack_corrected(:,:,i)=imwarp(ImageStack_corrected(:,:,i),TransformMatrix,'OutputView',imref2d(size(ImageStack_corrected(:,:,1))));            
                fTIF.WriteIMG(uint16(ImageStack_corrected(:,:,i)'+32768));
                disp(['Unwrapping and saving image number ',num2str(i)]);
            end 
            fTIF.close;
            disp(['Distortion correction finished, data saved: ',filename(1:end-4),'_wrappiing corrected.tif']);
            clear ImageStack_corrected
            
end