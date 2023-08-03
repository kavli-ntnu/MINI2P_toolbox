% Most Software Machine Data File

%% scanimage.SI (ScanImage)

% Global microscope properties
objectiveResolution = 75.082;     % Resolution of the objective in microns/degree of scan angle

% Data file location

% Custom Scripts
startUpScript = '';     % Name of script that is executed in workspace 'base' after scanimage initializes
shutDownScript = '';     % Name of script that is executed in workspace 'base' after scanimage exits

fieldCurvatureZs = [0 -200];     % Field curvature for mesoscope
fieldCurvatureRxs = [5000 5000];     % Field curvature for mesoscope
fieldCurvatureRys = [10000 10000];     % Field curvature for mesoscope
useJsonHeaderFormat = false;     % Use JSON format for TIFF file header

fieldCurvatureTip = 5;
fieldCurvatureTilt = 0;

%% scanimage.components.Motors (SI Motors)
% SI Stage/Motor Component.
motorXYZ = {'ScopeHolder' 'ScopeHolder' 'ScopeHolder'};     % Defines the motor for ScanImage axes X Y Z.
motorAxisXYZ = [2 1 3];     % Defines the motor axis used for Scanimage axes X Y Z.
scaleXYZ = [1 -1 1];     % Defines scaling factors for axes.
backlashCompensation = [0 0 0];     % Backlash compensation in um (positive or negative)

%% scanimage.components.scan2d.RggScan (MINI2P_001)

acquisitionDeviceId = 'vDAQ0';     % RDI Device ID
acquisitionEngineIdx = 1;

resonantScanner = 'MEMS2000Hz';     % Name of the resonant scanner
xGalvo = '';     % Name of the x galvo scanner
yGalvo = 'MEMS2000Hz_slowaxis';     % Name of the y galvo scanner
beams = {'920nm AOM'};     % beam device names
fastZs = {'TLens'};     % fastZ device names
shutters = {'Laser' 'MEMS' 'PMT'};     % shutter device names

channelsInvert = [true true false false];     % Logical: Specifies if the input signal is inverted (i.e., more negative for increased light signal)
keepResonantScannerOn = false;     % Always keep resonant scanner on to avoid drift and settling time issues

externalSampleClock = false;     % Logical: use external sample clock connected to the CLK IN terminal of the FlexRIO digitizer module
externalSampleClockRate = 8e+07;     % [Hz]: nominal frequency of the external sample clock connected to the CLK IN terminal (e.g. 80e6); actual rate is measured on FPGA
externalSampleClockMultiplier = 1;     % Multiplier to apply to external sample clock

extendedRggFov = 0;     % If true and x galvo is present, addressable FOV is combination of resonant FOV and x galvo FOV.

% Advanced/Optional
PeriodClockDebounceTime = 1e-07;     % [s] time the period clock has to be stable before a change is registered
TriggerDebounceTime = 5e-07;     % [s] time acquisition, stop and next trigger to be stable before a change is registered
reverseLineRead = 0;     % flips the image in the resonant scan axis
defaultFlybackTimePerFrame = 0.001;     % [s] default time to allow galvos to fly back after one frame is complete. overridden by cfg file
defaultFlytoTimePerScanfield = 0.001;     % [s] time to allow galvos to fly from one scanfield to the next. overridden by cfg file

% Aux Trigger Recording, Photon Counting, and I2C are mutually exclusive

% Aux Trigger Recording
auxTriggersTimeDebounce = 1e-07;     % [s] time after an edge where subsequent edges are ignored
auxTriggerLinesInvert = [false false false false];     % [logical] 1x4 vector specifying polarity of aux trigger inputs
auxTrigger1In = '';     % Digital input lines for aux trigger 1
auxTrigger2In = '';     % Digital input lines for aux trigger 2
auxTrigger3In = '';     % Digital input lines for aux trigger 3
auxTrigger4In = '';     % Digital input lines for aux trigger 4

% Signal Conditioning
disableMaskDivide = [false false false false];     % disable averaging of samples into pixels; instead accumulate samples

% I2C
i2cEnable = false;
i2cSdaPort = '';
i2cSclPort = '';
i2cAddress = 0;     % [byte] I2C address of the FPGA
i2cDebounce = 1e-07;     % [s] time the I2C signal has to be stable high before a change is registered
i2cStoreAsChar = false;     % if false, the I2C packet bytes are stored as a uint8 array. if true, the I2C packet bytes are stored as a string. Note: a Null byte in the packet terminates the string
i2cSendAck = true;     % When enabled FPGA confirms each packet with an ACK bit by actively pulling down the SDA line

% Laser Trigger
LaserTriggerPort = '';     % Digital input where laser trigger is connected.

% Trigger Outputs
frameClockOut = '/vDAQ0/D3.7';     % Output line for the frame clock
lineClockOut = '/vDAQ0/D0.4';     % Output line for the line clock
beamModifiedLineClockOut = '';     % Output line for beam clock
volumeTriggerOut = '';     % Output line for the volume clock

% Calibration data
scannerToRefTransform = [1 0 0;0 1 0;0 0 1];
virtualChannelSettings = [];
LaserTriggerDebounceTicks = 1;

enableHostPixelCorrection = false;
hostPixelCorrectionMultiplier = 500;

%% dabs.mirrorcle.ResonantAxis (MEMS2000Hz)
hAOZoom = '/vDAQ0/AO1';     % Analog output channel to be used to control resonant axis of the mirror e.g. '/vDAQ0/AO0'
hDOSync = '/vDAQ0/D1.7';     % String identifying the digital output for the sync signal e.g. '/vDAQ0/D0.1
hDOFilterX = '/vDAQ0/D3.2';     % String identifying the digital output for the X filter clock e.g. '/vDAQ0/D0.2
hDOFilterY = '/vDAQ0/D3.1';     % String identifying the digital output for the Y filter clock e.g. '/vDAQ0/D0.3

inputVoltageRange_Vpp = 3.87;     % Max input voltage range of the controller
angularRange_deg = 16;     % Max angular range of the device

% Default scan settings
syncPhase_deg = 120;
rampTime_s = 0;
xFilterClockFreq_Hz = 250000;
yFilterClockFreq_Hz = 100000;
xFilterClockEnable = 1;
yFilterClockEnable = 1;
nominalFrequency_Hz = 2000;

% Calibration Settings
amplitudeToLinePhaseMap = [0.8 8.4e-07;1.5 -2.664e-06;1.6 8.16e-06;2 -2e-06;3 -3.328e-06;4 -2.608e-06;4.5 -3.144e-06;4.615 -3.56e-06;5 -3.72e-06;5.333 -3.48e-06;5.455 -3.672e-06;5.714 -2.56e-06;6 -3.168e-06;6.154 -3.48e-06;6.667 -3.352e-06;7.273 -5.7912e-05;7.5 -3.416e-06;7.692 -3.504e-06;7.727 -3.416e-06;7.72727 0;8 -4.992e-06;8.182 6.3208e-05;8.333 -3.472e-06;8.5 -3.728e-06;8.636 -3.48e-06;9 -3.576e-06;9.091 -3.48e-06;9.5 -3.576e-06;10 0.00011808;11 -3.072e-06];     % translates an amplitude (degrees) to a line phase (seconds)
amplitudeLUT = zeros(0,2);     % translates a nominal amplitude (degrees) to an output amplitude (degrees)

%% dabs.generic.DigitalShutter (Laser)
DOControl = '/vDAQ0/D0.0';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = false;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% dabs.generic.DigitalShutter (MEMS)
DOControl = '/vDAQ0/D0.1';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = false;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% dabs.generic.DigitalShutter (PMT)
DOControl = '/vDAQ0/D0.2';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = false;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% dabs.generic.DigitalShutter (LED)
DOControl = '/vDAQ0/D0.3';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = false;     % invert output drive signal to shutter
openTime_s = 0.1;     % settling time for shutter in seconds

%% dabs.generic.BeamModulatorFastAnalog (920nm AOM)
AOControl = '/vDAQ0/AO2';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '/vDAQ0/AI0';     % feedback terminal e.g. '/vDAQ0/AI0'

outputRange_V = [0 2];     % Control output range in Volts
feedbackUsesRejectedLight = false;     % Indicates if photodiode is in rejected path of beams modulator.
calibrationOpenShutters = {'Laser'};     % List of shutters to open during the calibration. (e.g. {'Shutter1' 'Shutter2'}

powerFractionLimit = 1;     % Maximum allowed power fraction (between 0 and 1)

% Calibration data
powerFraction2ModulationVoltLut = [0 0;0.014 0.2;0.052 0.4;0.114 0.6;0.204 0.8;0.3013 1;0.4207 1.2;0.5538 1.4;0.6923 1.6;0.8461 1.8;1 2];
powerFraction2PowerWattLut = [0 0;1 0.169];
powerFraction2FeedbackVoltLut = zeros(0,2);
feedbackOffset_V = 0;

% Advanced Settings. Note: these settings are unused for vDAQ based systems
modifiedLineClockIn = '';     % Terminal to which external beam trigger is connected. Leave empty for automatic routing via PXI/RTSI bus
frameClockIn = '';     % Terminal to which external frame clock is connected. Leave empty for automatic routing via PXI/RTSI bus
referenceClockIn = '';     % Terminal to which external reference clock is connected. Leave empty for automatic routing via PXI/RTSI bus
referenceClockRate = 1e+07;     % if referenceClockIn is used, referenceClockRate defines the rate of the reference clock in Hz. Default: 10e6Hz

calibrationNumPoints = 40;
calibrationNumRepeats = 5;
calibrationAverageSamples = 5;
calibrationSettlingTime_s = 0.1;
calibrationFlybackTime_s = 0.1;

%% dabs.legacy.motor.LegacyMotor (ScopeHolder)
% Motor used for X/Y/Z motion, including stacks.

controllerType = 'PI C884/C863';     % If supplied, one of {'sutter.mp285', 'sutter.mpc200', 'thorlabs.mcm3000', 'thorlabs.mcm5000', 'scientifica', 'pi.e665', 'pi.e816', 'npoint.lc40x', 'bruker.MAMC'}.
comPort = '';     % Integer identifying COM port for controller, if using serial communication
customArgs = {};     % Additional arguments to stage controller. Some controller require a valid stageType be specified
invertDim = '$$$''---''';     % string with one character for each dimension specifying if the dimension should be inverted. '+' for normal, '-' for inverted
positionDeviceUnits = [-0.001 -0.001 -0.001];     % 1xN array specifying, in meters, raw units in which motor controller reports position. If unspecified, default positionDeviceUnits for stage/controller type presumed.
velocitySlow = [];     % Velocity to use for moves smaller than motorFastMotionThreshold value. If unspecified, default value used for controller. Specified in units appropriate to controller type.
velocityFast = [];     % Velocity to use for moves larger than motorFastMotionThreshold value. If unspecified, default value used for controller. Specified in units appropriate to controller type.
moveCompleteDelay = 0.1;     % Delay from when stage controller reports move is complete until move is actually considered complete. Allows settling time for motor
moveTimeout = '';     % Default: 2s. Fixed time to wait for motor to complete movement before throwing a timeout error
moveTimeoutFactor = '';     % (s/um) Time to add to timeout duration based on distance of motor move command

%% dabs.generic.FastZPureAnalog (TLens)
AOControl = '/vDAQ0/AO3';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '';     % feedback terminal e.g. '/vDAQ0/AI0'
FrameClockIn = '';     % frame clock input terminal e.g. '/Dev1/PFI0'

parkPositionUm = 0;     % park position in micron
travelRangeUm = [-250 10];     % travel range in micron

voltsPerUm = -0.032;     % volts per micron
voltsOffset = 0;     % volts that sets actuator to zero position

% Calibration Data
positionLUT = [-200 -190.9;-185 -176.68;-175 -167.2;-170 -160.12;-150 -142.2;-140 -135.84;-125 -123.7;-115 -116.3;-100 -102.6;-75 -84.1;-50 -65.6;-25 -43.2;-20 -35.34;0 0;10 26.77];     % Position LUT
feedbackVoltLUT = zeros(0,2);     % [Nx2] lut translating feedback volts into position volts

%% dabs.generic.GalvoPureAnalog (MEMS2000Hz_slowaxis)
AOControl = '/vDAQ0/AO0';     % control terminal  e.g. '/vDAQ0/AO0'
AOOffset = '';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '';     % feedback terminal e.g. '/vDAQ0/AI0'

angularRange = 8;     % total angular range in optical degrees (e.g. for a galvo with -20..+20 optical degrees, enter 40)
voltsPerOpticalDegrees = -1.1;     % volts per optical degrees for the control signal
parkPosition = 0;     % park position in optical degrees
slewRateLimit = Inf;     % Slew rate limit of the analog output in Volts per second

% Calibration settings
feedbackVoltLUT = zeros(0,2);     % [Nx2] lut translating feedback volts into position volts
offsetVoltScaling = 1;     % scalar factor for offset volts

voltsOffset = 0;

%% dabs.mini2p.MINI2P (MINI2P tools)
transformMatrixDirectory = 'C:\SI settings\2020a\TransformMatrix_mini2p2023_001.mat';     % transform Matrix Directory terminal  e.g. 'C:\Users\'
system = 'Forrest';     % Name of the system  e.g. 'System_001'
scope = 'MINI2P_2023_001';     % Name of the scope  e.g. 'MINI2P-L_001'
objective = 'D0309';     % Type of the objective  e.g. 'D0233'

