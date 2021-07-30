 %% Movie S8
 
aviobj1 = VideoWriter('MOVIE6.mp4','MPEG-4');
aviobj1.Quality=6;
aviobj1.FrameRate=50;
open(aviobj1)
%
Tiff_Files=dir(['D:\Science Movies\Movie S6\Final\*.tif']);
[Tiff_num ~]=size(Tiff_Files);
%%
for i=1:1:Tiff_num
    A=imread(['D:\Science Movies\Movie S6\Final\',Tiff_Files(i).name]);
    writeVideo(aviobj1,A);
    i
end
%
close (aviobj1)