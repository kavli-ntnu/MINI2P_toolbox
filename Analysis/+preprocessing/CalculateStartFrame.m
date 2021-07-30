function Startframe=CalculateStartFrame(FrameinEachSession, ImagingPlane)

FrameinEachSession=FrameinEachSession/ImagingPlane;

Startframe=zeros(1,size(FrameinEachSession,2));

if size(FrameinEachSession,2)==1
    Startframe=1;
else
    Startframe(1)=1;
    for i=2:1:size(FrameinEachSession,2)
        Startframe(i)=sum(FrameinEachSession(1:i-1))+1;
    end
end
end