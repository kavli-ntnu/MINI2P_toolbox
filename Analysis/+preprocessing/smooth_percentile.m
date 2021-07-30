function F_zero=smooth_percentile(F,TimeWindow,percentile)

Length=size(F,1);
F_zero=F;



HalfSpan=ceil(TimeWindow/2);

for i=1+HalfSpan:1:Length-HalfSpan
F_zero(i)=prctile(F(i-HalfSpan:1:i+HalfSpan),percentile);
end

for i=1:1:HalfSpan
F_zero(i)=prctile(F(i:1:i+2*HalfSpan),percentile);
end

for i=Length-HalfSpan+1:1:Length
F_zero(i)=prctile(F(i-2*HalfSpan:1:i),percentile);
end
end