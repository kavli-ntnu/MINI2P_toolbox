# MINI2P_toolbox 
 
## Project:

MINI2P freely-moving calcium imaging and spatial tuning analysis

## Description: 

MINI2P_toolbox includes the codes, softwares, 3D models, protocols, and etc. for buidling and using MINI2P to do freely-moving recording. MINI2P is an open-source miniature 2-photon microsocpe for fast high-resolution calcium imaging in freely-moving mice, published in Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". With the materials provided in this toolbox, people can assemble, test the MINI2P system, set up the animal tracking system, process the MINI2P imaging data, extract neuronal activity from single cells, and combine the neuronal activity data and the tracking data together for user-depedent downsteam analysis. The multi-FOV stitching software is also included. The codes for most of the anaylsis (grid cells, place cells, etc) in the paper "Large-scale two-photon calcium imaging in freely moving mice (2021)" are also provided.

## Contents: 

1) [Hardware](Hardware)

      [3D models (and 2D drawings for custom components)](Hardware) of all components for bulding a complete MINI2P system.
![image](https://user-images.githubusercontent.com/43905023/127703645-a6ea03ea-c1aa-4eaa-a9fd-1e6e75a082ed.png)

2) [Protocols](Protocols)

    a) [P1–Shopping and Machining list](https://github.com/kavli-ntnu/MINI2P_toolbox/blob/960fcab65f2849e5c8f486b81061c7c759f8e879/Protocols/P1%20-%20Shopping%20&%20Machining%20List%20.pdf). This document lists each essential component with the supplier, the product name, the model (or item reference) and its approximate price in Euro. 2D Drawings and 3D models of most components are available [here](https://github.com/kavli-ntnu/MINI2P_toolbox/tree/main/Hardware) .

    b) [P2-System building protocol](https://github.com/kavli-ntnu/MINI2P_toolbox/blob/main/Protocols/P2%20-%20System%20building%20protocol%20.pdf). This protocol includes all steps to assemble a MINI2P system. Each Protocol starts with a short-list of main reagents and tools needed, followed by an overview schematic of the module, and  a table with the main products.
       HC-920 assembly building video tutorial can be found on the link: TBD.
       Laser coupling video tutorial can be found on the link: TBD

    c) P3–Miniscope assembly protocol. How to assemble a MINI2P microscope is described in this protocol. Assembly tutorial video can be found on the [link](https://youtu.be/I0aYfi8GrIc).  Protocol document will be provided soon.

    d) P4–System operation manual. This section includes guidance for software installation, and imaging with MINI2P system. This protocol document will be provided soon.

    e) P5–System performance test protocol. This section illustrates how to calibrate and characterize the MINI2P system. This protocol document will be provided soon.

3) [Software](Software) 
      
      a) One scanimage Machine Data File (MDF) for [2000Hz MEMS-L scanner](Software/SI%20settings/Machine_Data_File_2000Hz.m). A second file for [5600Hz MEMS-F scanner](https://github.com/kavli-ntnu/MINI2P_toolbox/blob/main/Software/SI%20settings/Machine_Data_File_5600Hz.m)
      
      b) An example [Suite2P settings](Software/Suite2P%20options/GCaMP6S_P2_C1_7.25Hz_MEC.npy).
      
      c) Two DLC model configuration files:
      * [DLC1](Software/DLC%20model%20options/DLC1.yaml),
      * [DLC2](Software/DLC%20model%20options/DLC2.yaml),
      
      More details in [Documents](Documents/DeepLabCut-trained-Models.md).
      
      d) [AnimalTracker.vi](Software/AnimalTracker): a Labview program for recording animal behaviors and synchronizing the tracking camera recording with the MINI2P imaging. More details in [Documents](Documents/AnimalTracker.vi.md)

4) [Analysis](Analysis)

      a) [Pipelines](/Analysis/Pipeline) for spatial tuning analysis included in the paper (grid cells, place cells,etc).

      b) [NATEX.mlapp](/Analysis/Applications/NATEX): Nat Explorer, an application to load, process and preview the neuronal activity data (from the Suite2P output) and the tracking data (from the DLC output). It also combines the neuronal activity data and tracking data into the NAT.mat (Neuron Activity aligned with Tracking Matrix) and put all necessasy information into ExperimentInformation.mat for the user-specific downsteam analysis. More details in [Documents](Documents/NATEX.mlapp.md) 
      ![image](Analysis/Applications/NATEX/NATEX%20operation_speedup.gif)

      c) [StitchingChecker.mlapp](Analysis/Applications/StitchingChecker): an application to stitch multiple FOV recorded from different positions of the cortext. It can load in wide-field image as a reference for FOV alignment and can also take the retinotopic mapping result in for identifying different visual cortices. The precise alginment of FOVs is confirmed by 
      
      * overlapping of the landmarks between FOVs and the wide-field image, or between neighbouring FOVS; 
        
      * peak cross-correlation between FOVs and the wide-field image, or between neighbouring FOVs;
        
      * overlapping of the repeated cells in neighbouring FOVS. We also found this application can be used to register imagings recorded in multiple days. 
          
      More details in [Documents](DocumentsStitchingChecker.mlapp.md)
     ![image](Analysis/Applications/StitchingChecker/StitchingChecker%20operation_overview.gif)

     d) [DistortionCleaner.mlapp](Analysis/Applications/DistortionCleaner)： an application to elimiate the scanning distortion of MINI2P imaging, calibrate FOV and pixel size, and generate transform matrix. More details in [Documents](Documents/DistortionCleaner.mlapp.md)
     ![image](https://user-images.githubusercontent.com/43905023/127650948-b8ef7cc8-8c40-49b2-b374-dba90cc2844a.png)

5) [Documents](Documents)

    a) [Requirements](Documents/requirements.md): A list of the non-optical components necessary to build and use a MINI2P system, including licencsed software requirements. 
    
    b) [Bill of Materials](Documents/bill-of-materials.docx): A list of the the components (by source and part name) to build the optical system.
    
    c) [How-to](Documents/readme.md): A set of more detailed how-to documentation about how some components of the system work. 



   
## Usage:

  Applications NATEX, StitchingChecker and DistortionCleaner were written with Matlab app designer. In order to use these software, please press "open” in the home toolstrip of Matlab, select the software, wait until the app designer interface pops out, and then press "run". Some details in how to use each application is provided under [Documents](Documents/readme.md)


## Credits: 

This repositoriy is created and maintained by Simon Ball (simon.ball@ntnu.no) and Marina Rodrigues Jorge (marina.jorge@ntnu.no) in [Moser lab](https://www.ntnu.edu/kavli/moser-group#/view/about) at Kavli Institute for Systems Neuroscience, but got inputs and test from all authors in the paper Zong, et al.,"Large-scale two-photon calcium imaging in freely moving mice (2021)". Sections of the analysis code are based on the [Behavioural Neurology Toolbox](https://bitbucket.org/cnc-ntnu/bnt), (c) Vadim Frolov 2018.

MINI2P is a complete open-source project, we encourage people use, test, modify and further develop this toolbox. If you have any question or suggestion, or find any bugs in the codes, please contact us! If you use the code or data, please cite us!
