function out=MR_calculatrion(Ratemap,angle_range,distance_range)
% Ratemap=FullNeuronBehaviorDataSet.ECOFiringRateMap{2,2};
% angle_range=[0 2*pi];
% distance_range=[0,FullNeuronBehaviorDataSet.ECOMaxDistance];
%%
out=struct;
out.angle_range=angle_range;
out.distance_range=distance_range;
%%
[distance_num,angel_num]=size(Ratemap);
out.Angel=(linspace(out.angle_range(1),out.angle_range(2),angel_num))';
out.Distance=(linspace(out.distance_range(1),out.distance_range(2),distance_num))';
%%
out.MR=0;
for i=1:1:angel_num
    for j=1:1:distance_num
        out.MR=out.MR+Ratemap(j,i)*exp(1i*out.Angel(i));
    end
end
out.MR=out.MR./angel_num*distance_num;
%%
out.MRL=abs(out.MR);
out.MRA=angle(out.MR);
% MRA=MRA*180/pi;
if out.MRA<0
    out.MRA=2*pi+out.MRA;
else
end
%%
[~,out.Angle_nearest_ID]=min(abs(out.Angel-out.MRA));
out.Angle_nearest_Value=out.Angel(out.Angle_nearest_ID);
out.Rasponce_Distance_inPreferAngle=Ratemap(:,out.Angle_nearest_ID);
out.Distance_nearest_ID=find(out.Rasponce_Distance_inPreferAngle==max(out.Rasponce_Distance_inPreferAngle));
out.Distance_nearest_Value=out.Distance(out.Distance_nearest_ID);
%%
out.Distance_FWHM=out.Distance(preprocessing.FWHM(out.Rasponce_Distance_inPreferAngle));
%%
out.Rasponce_Angle_inPreferDistance=Ratemap(out.Distance_nearest_ID,:)';


end
