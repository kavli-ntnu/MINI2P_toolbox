% Delete cache that had been created for input file
%
function deleteCache(inputFile)
    global gBntData;
    try
        trials = io.parseGeneralFile(inputFile);
    catch
        io.readLinearInputFile(inputFile);
        trials = gBntData;
    end
    for i = 1:length(trials)
        trial = trials{i};

        [baseFolder, firstName] = helpers.fileparts(trial.sessions{1});
        if isempty(baseFolder)
            baseFolder = trial.path;
        end

        trial = io.detectRecordingSystem(trial);
        if strcmpi(trial.system, bntConstants.RecSystem.Neuralynx)
            baseFolder = trial.sessions{1};
        end
        
        if ~isfield(trial, 'basename') || isempty(trial.basename)
            trial.basename = firstName;
            trial.path = baseFolder;

            if length(trial.sessions) > 1
                if length(trial.sessions) > size(trial.cuts, 2)
                    [~, lastName] = helpers.fileparts(trial.sessions{end});
                    trial.basename = strcat(firstName, '+', lastName(end-1:end));
                end
            end
        end
        
        %% do the inheritance check
        % units
        if isempty(trial.units)
            if trial.inheritUnits
                if i > 1
                    trial.units = trials{i-1}.units;
                    trials{i}.units = trials{i-1}.units;
                end
            else
                if strcmpi(trial.system, bntConstants.RecSystem.Neuralynx)
                    trial = io.neuralynx.detectTetrodes(trial);
                else
                    trial = io.axona.detectTetrodes(trial);
                end
            end
        end
        if isempty(trial.units)
            error('BNT:noUnits', 'There are no units assigned with trial %u (first session ''%s'').\nCheck your input file!', i, trials{i}.sessions{1});
        end
        
        % shape info
        if ~isfield(trial.extraInfo, 'shape')
            if i > 1 && isfield(trials{i-1}.extraInfo, 'shape')
                trial.extraInfo.shape = trials{i-1}.extraInfo.shape;
                trials{i}.extraInfo.shape = trials{i-1}.extraInfo.shape;
            end
        end

        trial = io.initCuts(trial);

        tetrodes = unique(trial.units(:, 1), 'stable');
        numTetrodes = length(tetrodes);
        for t = 1:numTetrodes
            if t <= size(trial.cuts, 1)
                cutNames = trial.cuts(t, :);
            else
                if ~isempty(trial.cuts)
                    cutNames = trial.cuts(1, :);
                else
                    cutNames = {[]};
                end
            end
            notEmptyInd = ~cellfun(@isempty, cutNames);

            if ~isempty(trial.cuts) && all(notEmptyInd)
                unitCut = cutNames;
            else
                unitCut = [];
            end
            key = helpers.uniqueKey(trial.basename, unitCut);
            files = fullfile(trial.path, sprintf('*_T%u*_%s.mat', tetrodes(t), key));
            delete(files);
        end

        key = helpers.uniqueKey(trial.basename, trial.cuts);
        files = fullfile(trial.path, sprintf('*_%s.mat', key));
        delete(files);

        posFile = helpers.uniqueFilename(trial, 'posCleanScale');
        if exist(posFile, 'file')
            delete(posFile);
        end

        posFile = helpers.uniqueFilename(trial, 'posClean');
        if exist(posFile, 'file')
            delete(posFile);
        end

        posFile = helpers.uniqueFilename(trial, 'pos');
        if exist(posFile, 'file')
            delete(posFile);
        end
        % here is a small hack because the knowledge of kalmanFile is used
        [filePath, fileName] = helpers.fileparts(posFile);
        id = strfind(fileName, '_');
        fileName(id:end) = [];
        kalmanFile = fullfile(filePath, sprintf('%s_kalman*.mat', fileName));
        if exist(kalmanFile, 'file')
            delete(kalmanFile);
        end
    end
end