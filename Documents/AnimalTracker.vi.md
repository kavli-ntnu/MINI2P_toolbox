## Program location: 
MINI2P_toolbox/Software/AnimalTracker/
## Description: 
AnimalTracker.vi is a Labview program for recording animal behaviors, and synchronizing the tracking to the MINI2P imaging using a Basler camera. Currently, it was well validated by using camera acA2040-90um-NIR (https://www.baslerweb.com/en/products/cameras/area-scan-cameras/ace/aca2040-90um/). More validation will be published soon.
## Hardware requirement: 
1) Basler camera (recommend: [aca2040-90um](https://www.baslerweb.com/en/products/cameras/area-scan-cameras/ace/aca2040-90um/));
2) Basler camera [I/O cable](https://www.baslerweb.com/en/products/vision-components/cable/basler-gp-i-o-cable-hrs-6p-open-p-10-m/);
3) Basler cable [USB 3.0 Micro B](https://www.baslerweb.com/en/products/vision-components/cable/basler-cable-usb-3-0-micro-b-sl-a-p-3-m/);
## Software requirement: 
1) [NI-MAX](https://www.ni.com/en-no/support/downloads/drivers/download.system-configuration.html#371210)(20.0 or newer)
2) [Labview](https://www.ni.com/en-no/support/downloads/software-products/download.labview.html#369643)(2020 or newer)];
3) [NI-IMAQ](https://www.ni.com/en-no/support/downloads/drivers/download.vision-acquisition-software.html#367318);
3.1) [NI-IMAQmx](https://www.ni.com/en-no/support/downloads/drivers/download.ni-daqmx.html#409845);
3.2) [NI-IMAQdx](https://www.ni.com/en-no/support/downloads/drivers/download.vision-acquisition-software.html#367318);
4) (optional)[Pylon](https://www.baslerweb.com/en/sales-support/downloads/software-downloads/software-pylon-6-2-0-windows/) (6.2.0 or newer);
## Wiring:
To establish the correct camera wiring, first select Pins 2 and 5 of I/O cable which are Opto-coupled I/O input line and its ground, respectively. On the I/O cable, pins 2 and 5 are assigned as Pink and Grey, respectively. Link them with a BNC connector and plug the I/O cable from the camera to the vDAQ frame clock channel (one of the digital channels, in this case D3.7). See Fig. (1) for illustration.
![Schematic of I/O cable pins and labelling used to connect the Animal Tracker camera](https://github.com/WeijianZong/MINI2P_toolbox/blob/8b09d5eff4330c718e8ab2673bfa98cee6049cc3/Software/AnimalTracker/CAMERA-IOcable.png)

The Basler camera is working in External Trigger Mode, meaning it starts as soon as it gets a signal (0-5V) via input cable from the vDAQ (Start focus or grab in Scan Image, which sets up the Animal Tracker).
## Installation: 
The necessary drivers are those mentioned in Software Requirements from NI-MAX and NI-IMAQ. Pylon drivers are not needed, unless one wishes to visualize the camera using Pylon Software. In that case, USB 3.0 drivers should be installed.
### main pragram
To start-up, open the Virtual instrument application (.vi file) entitled "AnimalTracker.vi", and a graphical user interface will display as follows:


## GUI:
![Schematic of Animal Tracker GUI overview]
(https://github.com/kavli-ntnu/MINI2P_toolbox/blob/4a45033ba603d2f39251690c86131b548938de0d/Software/AnimalTracker/AnimalTracker-overview.png)

## Usage:
### change the default frame size, frame rate, exposure time and gain
1. On the same Directory where the imaging files are going to be saved, create an empty text file (in this case, "1.txt"). 
2. On the top left panel (see nr 1), select that text file created.
3. Chose Trigger Mode "External" and update frame rate (14.5 Hz) as well as other parameters as shown.
4. Finally, on the left bottom part, choose the camera type "1" for bottom and "2" for Top, and other relevant information which is going to be attached to a .csv file as:

###
## Output:

## .tiff to .mp4 conversion
