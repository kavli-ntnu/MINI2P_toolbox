%% MOVIE S9, read in vidoes
v = VideoReader('MOVIE_S9.mp4');
close all
TitleWords=zeros(720,1280,3);
TitleWords=insertText(TitleWords,[50,330],'Movie S9: example grid cells recorded from two planes in medial entorhinal cortex.','FontSize',28,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
% imshow(TitleWords)
Frame1=read(v,1);
aviobj1 = VideoWriter('MOVIE_S9_final.mp4','MPEG-4');
aviobj1.Quality=35;
aviobj1.FrameRate=75;
open(aviobj1)
SpendTime=4;
for i=1:1:SpendTime*75
    writeVideo(aviobj1,TitleWords/256);
    i
end
%
SpendTime=15;

for j=1:1:75*1
    Frame1=read(v,1);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=1:1:SpendTime*75
    Frame1=read(v,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end


for j=1:1:75*1
    Frame1=read(v,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x50','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    i
end

for i=SpendTime*75+1:10:16000
    Frame1=read(v,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end

for j=1:1:75*1
    Frame1=read(v,16001);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=16002:1:17000
    Frame1=read(v,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end
close (aviobj1)
%% MOVIE S10, read in vidoes
v2 = VideoReader('MOVIE_S10.mp4');
aviobj1 = VideoWriter('MOVIE_S10_final.mp4','MPEG-4');
aviobj1.Quality=42;
aviobj1.FrameRate=75;
open(aviobj1)
SpendTime=4;
TitleWords=zeros(600,1280,3);
TitleWords=insertText(TitleWords,[50,330],'Movie S10: example recording showing adjacent grid cells in medial entorhinal cortex.','FontSize',28,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);

for i=1:1:SpendTime*75
    writeVideo(aviobj1,TitleWords/256);
    i
end
%
SpendTime=15;

for j=1:1:75*1
    Frame1=read(v2,1);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=1:1:SpendTime*75
    Frame1=read(v2,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end


for j=1:1:75*1
    Frame1=read(v2,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x50','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);

    writeVideo(aviobj1,Frame1_corrected);
    i
end

for i=SpendTime*75+1:10:16000
    Frame1=read(v2,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end

for j=1:1:75*1
    Frame1=read(v2,16001);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=16002:1:17000
    Frame1=read(v2,i);
    Frame1_corrected=insertText(Frame1,[280,10],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end
close (aviobj1)

 %% MOVIE S11, read in vidoes
v4 = VideoReader('MOVIE_S11.mp4');
close all
TitleWords=zeros(720,1280,3);
TitleWords=insertText(TitleWords,[150,330],'Movie S11: example recording showing grid cells with tdTomato labeling','FontSize',28,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
TitleWords=insertText(TitleWords,[250,380],'       (projecting to DG/CA3) in medial entorhinal cortex.','FontSize',28,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
% imshow(TitleWords)
aviobj1 = VideoWriter('MOVIE_S11_final.mp4','MPEG-4');
aviobj1.Quality=37;
aviobj1.FrameRate=75;
open(aviobj1)
SpendTime=4;
for i=1:1:SpendTime*75
    writeVideo(aviobj1,TitleWords/256);
    i
end
%
SpendTime=15;

for j=1:1:75*1
    Frame1=read(v4,1);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
 	writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=1:1:SpendTime*75
    Frame1=read(v4,i);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
 	writeVideo(aviobj1,Frame1_corrected);
    i
end


for j=1:1:75*1
    Frame1=read(v4,i);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x50','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
 	   writeVideo(aviobj1,Frame1_corrected);
    i
end

for i=SpendTime*75+1:10:16000
    Frame1=read(v4,i);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x50','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);	
    writeVideo(aviobj1,Frame1_corrected);
    i
end


for j=1:1:75*1
    Frame1=read(v4,16001);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    Frame1_corrected=insertText(Frame1_corrected,[400,300],' Replay: time x5','FontSize',50,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',0.5);
    writeVideo(aviobj1,Frame1_corrected);
    j
end

for i=16002:1:17000
    Frame1=read(v4,i);
    Frame1_corrected=insertText(Frame1,[260,40],' Replay: time x5','FontSize',20,'TextColor',[256 256 256],'BoxColor',[0 0 0], 'BoxOpacity',1);
    writeVideo(aviobj1,Frame1_corrected);
    i
end



%
 close (aviobj1)
