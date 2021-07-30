function ImageStack_corrected=ImageUnwrapping(ImageStack,MovedPoint,RawPoint)


%% this function is for correcting the image distortion of the MINI2P images by using 2-D piecewise linear geometric transformation. The distortion is caused by the nonlinearity of the MEMS scanning.

%%
      ImageStack_corrected=ImageStack;
      tform = fitgeotrans(RawPoint,MovedPoint,'polynomial',3); %% generate the affine transformation mask
      Z_size=size(ImageStack,3);
      for j=1:1:Z_size
      ImageStack_corrected(:,:,j) = imwarp(ImageStack(:,:,j),tform,'OutputView',imref2d(size(ImageStack(:,:,j))));
%       ImageStack_corrected_16bit=uint16(ImageStack_corrected(:,:,j));
%       imwrite(ImageStack_corrected_16bit,[Result_folder,'\',FileName,'_wrappiing corrected.tif'],'WriteMode','append');
      disp(['Unwrapping image number ',num2str(j)]);
      end

end