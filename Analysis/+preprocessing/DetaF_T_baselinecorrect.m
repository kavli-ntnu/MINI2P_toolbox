% baseline correction based on the method reported in "Cellular resolution
% optical access to brain regions in fissures Imaging medial prefrontal
% cortex and grid cells in entorhinal cortex" PNAS,2014



function    [DetaF_F_filtered DetaF_raw F_zero Significant]=DetaF_T_baselinecorrect(F_raw,LowpassFrequency,MovingWindow,TimeThreadhold,Framerate,timeexpanded)

    DetaF_F_filtered=F_raw;
    DetaF_raw=F_raw;
    F_zero=F_raw; 
    F_zero_uncorrected=preprocessing.smooth_percentile(F_raw,ceil(MovingWindow*Framerate),8);
    F_std_local= movstd(F_raw,MovingWindow*Framerate);
    F_std_local_min=min(F_std_local);
    F_std_local_max=max(F_std_local);
    F_std_local_10percent=F_std_local_min+0.1*(F_std_local_max-F_std_local_min);
    F_lowpassed=lowpass(F_raw,LowpassFrequency,Framerate);
    F_M=mean(F_raw(find(F_std_local<F_std_local_10percent))-F_zero_uncorrected(find(F_std_local<F_std_local_10percent)));
    DetaF_F_filtered=(F_lowpassed-(F_zero_uncorrected+F_M))./(F_zero_uncorrected+F_M);
    DetaF_raw=(F_raw-(F_zero_uncorrected+F_M))./(F_zero_uncorrected+F_M);
    F_zero=F_zero_uncorrected+F_M;
    F_std_local_detaF=movstd(DetaF_F_filtered,TimeThreadhold*Framerate);    
    Significant_tem=(DetaF_F_filtered>2*F_std_local_detaF); %the significant threadhold was set as 1.5 X local STD
    Significant=zeros(size(DetaF_F_filtered));
    for j=1:1:length(Significant_tem)-ceil(TimeThreadhold*Framerate)-ceil(timeexpanded*Framerate)
        if sum(Significant_tem(j:1:j+ceil(TimeThreadhold*Framerate)))>=ceil(TimeThreadhold*Framerate)
            Significant(j:1:j+ceil(TimeThreadhold*Framerate)+ceil(timeexpanded*Framerate))=1;
        else 
        end
    end 
    for j=1+ceil(TimeThreadhold*Framerate)+ceil(timeexpanded*Framerate):1:length(Significant_tem)
        if sum(Significant_tem(j-ceil(TimeThreadhold*Framerate):1:j))>=ceil(TimeThreadhold*Framerate)
            Significant(j-ceil(TimeThreadhold*Framerate)-ceil(timeexpanded*Framerate):1:j)=1;
        else 
       end
    end
end