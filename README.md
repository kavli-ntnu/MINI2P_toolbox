# MINI2P_toolbox 
 
## Project:

MINI2P freely-moving calcium imaging and spatial tuning analysis

## Description: 

MINI2P_toolbox includes the codes, softwares, 3D models, protocols, and etc. for buidling and using MINI2P to do freely-moving recording. MINI2P is an open-source miniature 2-photon microsocpe for fast high-resolution calcium imaging in freely-moving mice, published in Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". With the materials provided in this toolbox, people can assemble, test the MINI2P system, set up the animal tracking system, process the MINI2P imaging data, extract neuronal activity from single cells, and combine the neuronal activity data and the tracking data together for user-depedent downsteam analysis. The multi-FOV stitching software is also included. The codes for most of the anaylsis (grid cells, place cells, etc) in the paper "Large-scale two-photon calcium imaging in freely moving mice (2021)" are also provided.

## Contents: 

1) [Hardware](https://github.com/kavli-ntnu/MINI2P_toolbox/tree/main/Hardware)

      [3D models (and 2D drawings for custom components)](Hardware/MINI2P_2021-07-12) of all components for bulding a complete MINI2P system.
![image](https://user-images.githubusercontent.com/43905023/127703645-a6ea03ea-c1aa-4eaa-a9fd-1e6e75a082ed.png)



2) [Software](Software) 
      
      a) Two scanimage Machine Data Files (MDF). One is for [2000Hz MEMS-L scanner](Software/SI%20settings/Machine_Data_File_2000Hz.m) and the other is for 5600Hz MEMS-F scanner;
      
      b) An example [Suite2P settings](https://github.com/WeijianZong/MINI2P_toolbox/blob/main/Software/Suite2P%20options/GCaMP6S_P2_C1_7.25Hz_MEC.npy).
      
      c) Three DLC model configuration files: [
      * DLC1](Software/DLC%20model%20options/DLC1.yaml),
      * [DLC2](Software/DLC%20model%20options/DLC2.yaml),
      * [DLC3](Software/DLC%20model%20options/DLC3.yaml)
      More details in [Wiki page](https://github.com/WeijianZong/MINI2P_toolbox/wiki/DeepLabCut-trained-Models).
      
      d) [AnimalTracker.vi](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Software/AnimalTracker): a Labview pragram for recording animal behaviors and synchronizing the tracking camera recording with the MINI2P imaging. More details in [Wiki page](https://github.com/WeijianZong/MINI2P_toolbox/wiki/AnimalTracker.vi)

3) [Analysis](Analysis)

      a) [Pipelines](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Analysis/Pipeline) for spatial tuning analysis included in the paper (grid cells, place cells,etc).

      b) [NATEX.mlapp](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Analysis/Applications/NATEX): Nat Explorer, an application to load, process and preview the neuronal activity data (from the Suite2P output) and the tracking data (from the DLC output). It also combines the neuronal activity data and tracking data into the NAT.mat (Neuron Activity aligned with Tracking Matrix) and put all necessasy information into ExperimentInformation.mat for the user-specific downsteam analysis. More details in [Wiki page](https://github.com/WeijianZong/MINI2P_toolbox/wiki/NATEX.mlapp) 
      ![image](https://github.com/WeijianZong/MINI2P_toolbox/blob/acecff10421ba8fe03a35a9fc02a8ed57fc91d2a/Analysis/Applications/NATEX/NATEX%20operation_speedup.gif)

      c) [StitchingChecker.mlapp](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Analysis/Applications/StitchingChecker): an application to stitch multiple FOV recorded from different positions of the cortext. It can load in wide-field image as a reference for FOV alignment and can also take the retinotopic mapping result in for identifying different visual cortices. The precise alginment of FOVs is confirmed by 
      
      * overlapping of the landmarks between FOVs and the wide-field image, or between neighbouring FOVS; 
        
      * peak cross-correlation between FOVs and the wide-field image, or between neighbouring FOVs;
        
      * overlapping of the repeated cells in neighbouring FOVS. We also found this application can be used to register imagings recorded in multiple days. 
          
      More details in [Wiki page](https://github.com/WeijianZong/MINI2P_toolbox/wiki/StitchingChecker.mlapp)
     ![image](https://github.com/WeijianZong/MINI2P_toolbox/blob/effa456578ddc0eb99795abb82c57117ab69400c/Analysis/Applications/StitchingChecker/StitchingChecker%20operation_overview.gif)

     d) [DistortionCleaner.mlapp](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Analysis/Applications/DistortionCleaner)：an application to elimiate the scanning distortion of MINI2P imaging, callibrate FOV and pixel size, and generate transform matrix. More details in [Wiki page](https://github.com/WeijianZong/MINI2P_toolbox/wiki/DistortionCleaner.mlapp)
     ![image](https://user-images.githubusercontent.com/43905023/127650948-b8ef7cc8-8c40-49b2-b374-dba90cc2844a.png)

4) [Documents](https://github.com/WeijianZong/MINI2P_toolbox/tree/main/Documents)


## Usage:

  Applications NATEX, StitchingChecker and DistortionCleaner were written with Matlab app designer. In order to use these softwares, please press "open” in the home toolstrip of Matlab, select the software, wait until the app designer interface pops out, and then press "run". The details about how to use these codes and softwares are in wiki pages of this repository.


## Credits: 

The codes were mainly written by Weijian Zong (weijian.zong@ntnu.no) in [Moser lab](https://www.ntnu.edu/kavli/moser-group#/view/about) at Kavli Institute for Systems Neuroscience, but got inputs and test from all authors in the paper Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". Sections of the analysis code are based on the [Behavioural Neurology Toolbox], (c) Vadim Frolov 2018.

MINI2P is a complete open-source project, we encourage people use, test, modify and further develop this toolbox. If you have any question or suggestion, or find any bugs in the codes, please contact us! If you use the code or data, please cite us!

