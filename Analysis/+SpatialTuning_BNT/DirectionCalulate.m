function hd = DirectionCalulate(X1,Y1,X2,Y2)
    hd=rem(atan2d(Y2-Y1,X2-X1)+360, 360);
end