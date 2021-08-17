function OverlapRatio= overlap(ROI1,ROI2,ImageSize)

Image=zeros(ImageSize);

Xpoint=mod(double(ROI1(1,:))',ImageSize(1));
Ypoint=mod(double(ROI1(2,:))',ImageSize(2));
Xpoint(Xpoint==0)=ImageSize(1);
Ypoint(Ypoint==0)=ImageSize(2);
for m=1:1:length(Xpoint)
      Image(Ypoint(m),Xpoint(m))=1;
end 

Xpoint=mod(double(ROI2(1,:))',ImageSize(1));
Ypoint=mod(double(ROI2(2,:))',ImageSize(2));
Xpoint(Xpoint==0)=ImageSize(1);
Ypoint(Ypoint==0)=ImageSize(2);
for m=1:1:length(Xpoint)
      Image(Ypoint(m),Xpoint(m))=Image(Ypoint(m),Xpoint(m))+1;
end 


OverlapRatio(1)=length(find(Image(:)==2))./size(ROI1,2);
OverlapRatio(2)=length(find(Image(:)==2))./size(ROI2,2);
OverlapRatio(3)=length(find(Image(:)==2))./(length(find(Image(:)==2))+length(find(Image(:)==1)));

end