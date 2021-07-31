% Initialize trial with default values
%
function trial = initTrial()
    % Recording system. Possible values are 'Axona' and 'NeuraLynx'.
    trial.system = nan;

    % Session base name. I.e. for single session '2103201101'. For combined sessions
    % it will be '2103201101+02'. More than two combined sessions are not supported.
    trial.basename = '';

    % Session path, normally the folder of the first session
    trial.path = '';

    trial.positions = []; % Matrix Nx5, N - number of entries.
                                % columns:
                                % 1 timestamps
                                % 2 LED1 x values
                                % 3 LED1 y values
                                % 4 LED2 x values, also known as x2
                                % 5 LED2 y values, also knows as y2

    trial.spikes = [];    % Matrix Nx3, N - number of entries. 1 column stores timestamps, spike data.
                                % 2 column contains tetrode number (integer), 3 column - unit
                                % number (integer).

    % Extra information, structure. Could contain fields
    % room - string with information about recording room;
    % shape - struct with information about recording box shape, fields:
    %         descr - raw string from file
    %         type - integer. 1 = box, 2 = cylinder, 3 = linear track
    %         value - double. Side length or diameter of the arena. If dimensions are not equal,
    %                 then value will have several values. For example for box arena value(1)
    %                 defines x-dimension, and value(2) defines y-dimension.
    % startTime - string
    % calibration - structure with following fields:
    %         file - Path to file with camera calibration information
    %         data - calibration data
    %
    % innerSize - structure with following fields:
    %         rec - Path to recording with LEDs that define arena inner size
    %         data - extracted inner size data
    % outerSize - structure with following fields
    %         rec - Path to recording with LEDs that define arena outer size
    %         data - extracted outer size data
    % subjectId - string, animal/subject ID
    %
    trial.extraInfo = struct();

    % Sampling time. Default values are from Axona system.
    trial.sampleTime = 0.02; % 50 Hz
    trial.videoSamplingRate = 50;

    trial.sessions = {}; % session names
    trial.cuts = {}; % cut filenames. Matrix of size NxM. Each row contains cut filenames
                     % for a particular tetrode (tetrodes are obtained from trial.units).
                     % Values in columns allow to load multiple cut files for cells,
                     % or have multiple cuts per tetrode in case of combined sessions.
    % Information about cut files that have been used to define units.
    % The key is an integer that represents the cell: each tetrode and cell number is
    % mapped in a single value by function helpers.tetrodeCellLinearIndex.
    % Each value is a Nx3 cell array. N is the number of cut files that have been used
    % to define this cell. First column contains path to cut file, second column
    % contains corresponding MD5 hash value. Both columns are strings.
    % Third column contains session index of a session to which this cut file belongs.
    % This is needed to adjust spike offset in case of multiple session loading.
    % Used to determine if cut files have been changed since last load of data in Matlab.
    trial.cutHash = containers.Map('KeyType', 'int64', 'ValueType', 'any');

    trial.tetrode = [];
    trial.units = []; % Matrix of units, Nx2. First column contains
                      % tetrode numbers, second - cell numbers.

    trial.startIndices = []; % When multiple sessions are combined, it can be
                             % useful to track position of each session individually.
                             % This vector stores starting index of each session. Index
                             % is applied to position matrix.

     % List of LFP channels to load/use. -1 has a special meaning being
     % 'all available channels'.
     trial.lfpChannels = [];
     
     % Flag that determines if this trial inherits units of a previous
     % trial or not. Default behaviour is that units are not inherited.
     trial.inheritUnits = false;
     
     % Identifier of cutting software. See also bntConstants.CutApp.
     trial.cutApp = [];
     
     % Fields that are going to be initialized later if at all:
     % lfpInfo - instance of io.LfpInfoInterface class
end
