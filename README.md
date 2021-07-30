# MINI2P_toolbox
 
## Project:
  MINI2P_toolbox

## Description: 
  MINI2P_toolbox includes the codes, softwares, 3D models, protocols, and etc. for buidling and using MINI2P to do freely-moving recording. MINI2P is an open-source miniature 2-photon microsocpe for fast high-resolution calcium imaging in freely-moving mice, published in Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". With the materials provided in this toolbox, people can assemble, test the MINI2P system, set up the animal tracking system, process the MINI2P imaging data, extract neuronal activity from single cells, and combine the neuronal activity data and the tracking data together for user-depedent downsteam analysis. The multi-FOV stitching software is also included. The codes for most of the anaylsis (grid cells, place cells, etc) in the paper "Large-scale two-photon calcium imaging in freely moving mice (2021)" are also provided.

## Contents: 
1) Hardware: 3D models (and 2D drawings for custom components) of all componets for bulding a complete MINI2P system.

2) Software: include a list of requiered software to run complete freely-moving MINI2P recordings, an example suite2P setting option, three DLC models (model 1 was trained by vidoes recorded from bottom of transparent floors with 'nose', 'mouse', 'lefthand', 'righthand', 'leftleg', 'rightleg', 'tailbase', and 'bodycenter' extracted; model 2 was trained by vidoes recorded from top of boxes with 'miniscope', 'leftear', 'rightear', 'bodycenter', and 'tailbase' extracted; model 3 was trained by vidoes recorded from top of boxes with 'nose', 'leftear', 'rightear', 'bodycenter', and 'tailbase' extracted), and a Labview-based traking vedio recording software.

3) Analysis: include pipelines for spatial tuning analysis included in the paper (grid cells, place cells,etc), and three Matlab-based applications: NATEX, StitchingChecker and DistortionCleaner.

      a) NATEX: Nat Explorer, an application to load, process and preview the neuronal activity data (from the Suite2P output) and the tracking data (from the DLC output). It                 also combines the neuronal activity data and tracking data into the NAT.mat (Neuron Activity aligned with Tracking Matrix) and put all necessasy information                     into ExperimentInformation.mat for the user-specific downsteam analysis.
      ![image](https://user-images.githubusercontent.com/43905023/127646675-27e2dd97-a133-4289-85f3-d8024797a86c.png)

      b) StitchingChecker: an application to stitch multiple FOV recorded from different positions of the cortext. It can load in wide-field image as a reference for FOV                                  alignment and can also take the retinotopic mapping result in for identifying different visual cortices. The precise alginment of FOVs is confirmed                              by  i) overlapping of the landmarks between FOVs and the wide-field image, or between neighbouring FOVS; ii) peak cross-correlation between FOVs and                            the wide-field image, or between neighbouring FOVs; and iii) overlapping of the repeated cells in neighbouring FOVS. We also found this application is                            can be used to register imagings recorded in multiple days.
     ![image](https://user-images.githubusercontent.com/43905023/127649125-7b5b9de3-53f1-4d65-8fc2-63f9bdea75a4.png)

     c) DistortionCleaner：an application to elimiate the scanning distortion of MINI2P imaging, callibrate FOV and pixel size, and generate transform matrix (see original paper for details).
     ![image](https://user-images.githubusercontent.com/43905023/127650948-b8ef7cc8-8c40-49b2-b374-dba90cc2844a.png)



## Installation: 
  Most of the codes in this toolbox were written in Matlab (the python version will be published in the future). Please make sure the Matlab (recommend version 2020b) was installed in your computer. Firstly, download or clone the whole repository in a local foler. Then, include folder and subfolders \MINI2P_toolbox\MINI2P_toolbox\codes in Matlab search path by pressing "set path” in the home toolstrip of Matlab, then “add with subfolders...". 

## Usage:
  Applications NATEX, StitchingChecker and DistortionCleaner were written with Matlab app designer. In order to use these softwares, please press "open” in the home toolstrip of Matlab, select the software, wait until the app designer interface pops out, and then press "run". The details about how to use these codes and softwares are in wiki pages of this repository.


## Credits: the codes were mainly written by Weijian Zong (weijian.zong@ntnu.no) at Kavli Institute for Systems Neuroscience, but got inputs and test from all authors in the paper Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". Since MINI2P is a complete open-source project, we encourage people use, test, modify and further develop this toolbox. If you have any question or suggestion, or find any bugs in the codes, please contact us! If the codes and software in this toolbox help your project, please cite our paper Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)" in your publications! Enjoy it!

