function WJ_DistortionCorrection_3D_SItriggered(src,evt,varargin)
            global filename        
            filename=[src.hSI.hScan2D.logFilePath,'\', src.hSI.hScan2D.logFileStem,'_',num2str(src.hSI.hScan2D.logFileCounter-1,'%05d'),'.tif'];                                   
            disp(['Just saved raw data: ',filename]);
            disp('Start the distortion correction....');
            disp('Reloading the raw data....');
     
%             filename='97045_20210308_ML-400_AL-400_1Openfiled_00001.tif';
            
            [header, RawImage, imgInfo] = scanimage.util.opentif(filename);
            filename1=imgInfo.filename;
            disp(['Rawe 2P imaging data: "' filename1  '" is loaded']);
            save([filename(1:end-4),'_tifHeader.mat'],'header');% save the SI information, only necessasy for SI information saving. No need for DJ analysis
            disp([filename(1:end-4),'_tifHeader.mat',' is saved']);
            save([filename(1:end-4),'_imgInfo.mat'],'imgInfo');% save the imaging information,only necessasy for SI information saving. No need for DJ analysis
            disp([filename(1:end-4),'_imgInfo.mat',' is saved']);
            %
            Dimention=size(RawImage); 
            Zs=header.SI.hStackManager.zs;
            % load in 3D transform matrx according to different frame szie
            if Dimention(1)==512
                % read 515 3D Tranform matrix
                TransformMatrix_tem=load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\TransMatrix_512_3D.mat');
                % read 512 3D FOV information
                FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\FOVreport_512_20211111.csv');
                Depth=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\Depth20211111.csv');
                %
            elseif Dimention(1)==256
                % read 256 3D Tranform matrix
                TransformMatrix_tem=load('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\TransMatrix_256_3D.mat');
                % read 256 3D FOV information
                FOV=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\FOVreport_256_20211111.csv');
                Depth=csvread('\\forskning.it.ntnu.no\ntnu\mh-kin\moser\open2pmini\FOV callibration\RubenWeijianNienke\MINI2P_L_WJ001\Report\Depth20211111.csv');

                %
            else
            end
            %
            TransformMatrix=struct2cell(TransformMatrix_tem);
            
            PlaneMark=[];
            if length(Dimention)<=3                  
                ImageStack_corrected=RawImage;
                PlaneMark=ones(size(ImageStack_corrected(:)),1);
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
                        PlaneMark(k)=1;
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
                                PlaneMark(k)=m;
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
                                PlaneMark(k)=m;
                                k=k+1;
                            end
                        end
                    end
                end    
            end
            clear RawImage
            % find which matrix to use
            PlaneClosest=PlaneMark;
            for i=1:1:size(PlaneMark,2)
                PlaneClosest(i)=Zs(PlaneMark(i));
                [minValue,closestIndex]=min(abs(Depth-PlaneClosest(i)));
                PlaneClosest(i)= closestIndex;
            end

            fTIF = DataIO.Fast_BigTiff_Write([filename(1:end-4),'_wrappiing-corrected.tif'],1*10000/FOV(3),0);
            for i=1:1:size(ImageStack_corrected,3)               
                ImageStack_corrected(:,:,i)=imwarp(ImageStack_corrected(:,:,i),TransformMatrix{PlaneClosest(i),1},'OutputView',imref2d(size(ImageStack_corrected(:,:,1))));            
                fTIF.WriteIMG(int16(ImageStack_corrected(:,:,i)'));
                disp(['Unwrapping and saving image number ',num2str(i)]);
            end 
            fTIF.close;
            disp(['Distortion correction finished, data saved: ',filename(1:end-4),'_wrappiing corrected.tif']);
%             clear ImageStack_corrected            
end