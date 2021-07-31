% uniqueFilename - Provide a file name used to store data on hard disc.
%
% Meta iformation and data are stored on hard disc in a file. In order not to mix
% up the data it is stored in a uniqely-named files. This function provides a name
% for different kind of files.
%
%  USAGE
%    name = uniqueFilename(sessionData, fileKind, ...);
%    sessionData        Structure with information about the current session. Important
%                       are fields .basename and .cuts.
%    fileKind           String. Determines the ouput name. Supported file kinds are:
%       cut             File with spike timestamps. Requires additional argument 'unit'. Unit is
%                       a two-element vector containing tetrode and cell numbers.
%                       Returned filename depends on session names and cut name.
%       posCleanScale   File with position data (cleaned and scaled). Not in use from 04.11.2014!
%       posClean        File with position data (cleaned).
%       data            File with information about a trial.
%       pos             File with raw position data. Filename depends only
%                       on original names of sessions, i.e. trial basename.
%
%    Additional args    Additional arguments depend on fileKind. This can be a unit number
%                       if fileKind is 'cut'.
%    name               Generated name.
%
function name = uniqueFilename(sessionData, fileKind, varargin)
    key = helpers.uniqueKey(sessionData.basename, sessionData.cuts);

    switch fileKind
        case 'cut'
            if nargin < 3
                error(message('MATLAB:minrhs'));
            end
            unit = varargin{1};

            if length(unit) == 1 && ~isfield(sessionData, 'tetrode')
                error('BNT:arg', 'Invalid additional argument for type ''cut''. The value should be two-elements vector, single value was provided.');
            end

            tetrodes = unique(sessionData.units(:, 1), 'stable');
            ind = find(tetrodes == unit(1));

            % in case user specified only one cut file, which is actually a pattern
            if ~isempty(ind)
                if ind <= size(sessionData.cuts, 1)
                    cutNames = sessionData.cuts(ind(1), :);
                else
                    if ~isempty(sessionData.cuts)
                        cutNames = sessionData.cuts(1, :);
                    else
                        cutNames = {[]};
                    end
                end
            else
                cutNames = {[]};
            end
            notEmptyInd = ~cellfun(@isempty, cutNames);

            if ~isempty(sessionData.cuts) && all(notEmptyInd)
                unitCut = cutNames;
            else
                unitCut = [];
            end
            key = helpers.uniqueKey(sessionData.basename, unitCut);

            if length(unit) > 1
                name = fullfile(sessionData.path, sprintf('%s_T%uC%u_%s.mat', ...
                    sessionData.basename, ...
                    unit(1), unit(2), key));
            else
                name = fullfile(sessionData.path, sprintf('%s_T%uC%u_%s.mat', ...
                    sessionData.basename, ...
                    sessionData.tetrode, unit, key));
            end
        case 'posCleanScale'
            name = fullfile(sessionData.path, sprintf('%s_posCleanScale.mat', ...
                sessionData.basename));

        case 'posClean'
            name = fullfile(sessionData.path, sprintf('%s_posClean.mat', ...
                sessionData.basename));

        case 'data'
            name = fullfile(sessionData.path, sprintf('%s_data_%s.mat', ...
                sessionData.basename, key));

        case 'pos'
            name = fullfile(sessionData.path, sprintf('%s_pos.mat', ...
                sessionData.basename));

        otherwise
            error('BNT:arg', 'Unknown value of ''fileKind'' argument');
    end
end
