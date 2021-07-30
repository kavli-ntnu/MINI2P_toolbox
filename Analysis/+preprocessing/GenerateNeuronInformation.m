function [NeuronInformation, NeuronActiveMatrix] =GenerateNeuronInformation(Suite2P_file)

Suite2P_raw=load(Suite2P_file);
disp(['Suite2P result: ',Suite2P_file,' is loaded.']);

NeuronInformation=struct;
NeuronInformation.TotalROINumber=length(Suite2P_raw.stat(1,:));
NeuronInformation.FrameClockNumber=length(Suite2P_raw.F(1,:));
NeuronInformation.ROIStat=Suite2P_raw.stat;
NeuronInformation.CellPossibility=Suite2P_raw.iscell;
% NeuronInformation.TotalIsCellROINumber=length(Suite2P_raw.iscell(:,2)==1);
NeuronInformation.Options=Suite2P_raw.ops;
NeuronInformation.volumerate=NeuronInformation.Options.fs; %Hz
NeuronInformation.ImagingPlane=NeuronInformation.Options.nplanes;
NeuronInformation.ImagingChannel=NeuronInformation.Options.nchannels;

NeuronActiveMatrix=struct;
NeuronActiveMatrix.EventTrain=zeros(NeuronInformation.FrameClockNumber,2+NeuronInformation.TotalROINumber);
NeuronActiveMatrix.EventTrain(:,3:end)=Suite2P_raw.spks';
NeuronActiveMatrix.F_raw=zeros(NeuronInformation.FrameClockNumber,2+NeuronInformation.TotalROINumber);
NeuronActiveMatrix.F_raw(:,3:end)=Suite2P_raw.F';
NeuronActiveMatrix.F_neuropil=zeros(NeuronInformation.FrameClockNumber,2+NeuronInformation.TotalROINumber);
NeuronActiveMatrix.F_neuropil(:,3:end)=Suite2P_raw.Fneu';
NeuronActiveMatrix.EventTrain(:,1)=0:1/NeuronInformation.volumerate:(length(NeuronActiveMatrix.EventTrain(:,1))-1)/NeuronInformation.volumerate;
NeuronActiveMatrix.F_raw(:,1)=0:1/NeuronInformation.volumerate:(length(NeuronActiveMatrix.EventTrain(:,1))-1)/NeuronInformation.volumerate;
NeuronActiveMatrix.F_neuropil(:,1)=0:1/NeuronInformation.volumerate:(length(NeuronActiveMatrix.EventTrain(:,1))-1)/NeuronInformation.volumerate;

disp(['NeuronActiveMatrix for: ',Suite2P_file,' was generated.']);
   
end