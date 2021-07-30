function Speed=SpeedCalculate(Position,timeinterval,smoothwidth)

    Position(:,1)=general.smoothGauss(Position(:,1),smoothwidth);
    Position(:,2)=general.smoothGauss(Position(:,2),smoothwidth);

    Speed=zeros(size(Position,1),1);

for i=1:1:size(Position,1)-1
Speed(i)=sqrt((Position(i+1,1)-Position(i,1))^2+(Position(i+1,2)-Position(i,2))^2)/timeinterval;
end
Speed(end)=Speed(end-1);

Speed=general.smoothGauss(Speed,smoothwidth);

end