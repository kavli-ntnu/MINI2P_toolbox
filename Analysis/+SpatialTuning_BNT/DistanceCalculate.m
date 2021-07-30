function Ds=DistanceCalculate(Position1,Position2)
Tem=Position2-Position1;
Ds=zeros(size(Tem,1),1);
for i=1:1:size(Ds)
    Ds(i)=sqrt(Tem(i,1).^2+Tem(i,2).^2);
end
end