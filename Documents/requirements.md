# Requirements for the mini2p toolbox

## Hardware

Substituting alternative, comparable, pieces of hardware is generally possible. Modifications to the software will usually be required if so.

### Microscope

A complete bill of materials can be provided upon request. Details of specific hardware components of the microscope are discussed in the associated publication:

> W. Zong *et al*, "Large-scale two-photon calcium imaging in freely moving mice". [BioRXiv (2021)](https://doi.org/10.1101/2021.09.20.461015)

A single computer workstation is required for experiment control. Experimental control is not particular demanding of computing resources:
* Windows 10 x64
* Spare PCIe 3 x8 interface
* Spare PCIe 6-pin power connector

Post-recording imaging analysis via Suite2P and DeepLabCut is more demanding. Post-recording analysis does not need to take place on the same computer. 
* Version 11 Cuda compatible graphics card with >8GB vRAM
* Minimum 4GB RAM per CPU core, recommended 8GB.
* Minimum 4 cores. 



## Software

### Recording

* Matlab 2020b or more recent
  * [Mathworks](https://se.mathworks.com/)
* LabView 2020 or more recent
  * [National Instruments](https://www.ni.com/en-no/shop/labview.html)
  * [Anaconda](https://www.anaconda.com/products/individual)
* ScanImage 2020 Basic or more recent
  * The Free edition does not provide the necessary tools
  * [Vidrio](https://vidriotechnologies.com/compare-scanimage-version/)
* vDAQ RDI driver 2019-11-20 or more recent
  * SDK for the vDAQ interface card
  * [Vidrio](https://vidriotechnologies.com/drivers)
* Pylon 6.2.0 or more recent
  * Camera controller software for the Basler camera series
  * [Basler](https://www.baslerweb.com/en/sales-support/downloads/software-downloads/)
* PMT2100 4.0 or more recent
  * SDK for the Thorlabs PMT2100 photomultiplier tube
  * [Thorlabs](https://www.thorlabs.com/software_pages/ViewSoftwarePage.cfm?Code=PMT)

### Postprocessing

* Matlab 2020b or more recent
  * [Mathworks](https://se.mathworks.com/)
* Python 3.8 or more recent
  * Due to Suite2P's requirements, Anaconda and the Conda package manager are recommended
* DeepLabCut 2.2.0 or more recent
  * MachineMarkerless subject tracking
  * [Mathis Lab](http://www.mackenziemathislab.org/deeplabcut)
* Suite2p v0.10.0 or more recent
  * Two-photon imaging processing
  * [Pachitaru Lab](https://github.com/MouseLand/suite2p)
