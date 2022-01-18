%% load in anchor pairs in different forcal plane (callibrated by 512x512)

AnchorPairlist=dir(['*.csv']);

AnchorPair=cell(size(AnchorPairlist,1),1);

for i=1:1:size(AnchorPairlist,1)
    AnchorPair{i,1}=csvread(AnchorPairlist(i,1).name);
end
%% knowledge from recording

Depth=[0:-10:-240]'; % Tlens position, um

TransMatrix_512_3D=cell(size(Depth,1),1);

for i=1:1:size(Depth,1)
    TransMatrix_512_3D{i,1}=fitgeotrans(AnchorPair{i,1}(:,3:4),AnchorPair{i,1}(:,1:2),'pwl');
end

%% test
select=1;
P1=imread('FOVStack00.tif');
P1=double(P1);
P1_corrected=imwarp(P1,TransMatrix_512_3D{select,1},'OutputView',imref2d(size(P1)));

close all

figure

subplot(1,2,1)

imshow(P1/max(P1(:)));

subplot(1,2,2)

imshow(P1_corrected/max(P1_corrected(:)));

%% anchor pairts of 512 to 256

AnchorPair_scaled=cell(size(AnchorPairlist,1),1);

for i=1:1:size(AnchorPairlist,1)
    AnchorPair_scaled{i,1}=AnchorPair{i,1}/2;
end
%% transfer matrix of 512 to 256

TransMatrix_256_3D=cell(size(Depth,1),1);

for i=1:1:size(Depth,1)
    TransMatrix_256_3D{i,1}=fitgeotrans(AnchorPair_scaled{i,1}(:,3:4),AnchorPair_scaled{i,1}(:,1:2),'pwl');
end
%% test the 256 image

select=1;
P1=imread('FOVStack19.tif');
P1=double(P1);
P1_corrected=[];
P1_corrected=imwarp(P1,TransMatrix_256_3D{select,1},'OutputView',imref2d(size(P1)));

close all

figure

subplot(1,2,1)

imshow(P1/max(P1(:)));

subplot(1,2,2)

imshow(P1_corrected/max(P1_corrected(:)));
%%
