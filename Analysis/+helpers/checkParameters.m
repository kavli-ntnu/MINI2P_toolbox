% Check that parameters are up to date
%
% The most recent parameters are found in file 'examples/bntSettings.m'. However,
% users are encouraged to have their own version of this file. This function
% helps to track down the difference between users's parameter and those of BNT.
%
%  USAGE
%   params = helpers.checkParameters(userParams)
%   userParams      Structure with parameters
%   params          If user parameters miss some values, then default
%                   values will be added. Combined set of parameters is
%                   returned.
%                   If user parameters set doesn't differ from default,
%                   then params = userParams.
%
%  OUTPUT
%   In case parameters are the same you will see no output.
%   Otherwise, the will be a warning message and suggestion to fix your parameters.
%
function p = checkParameters(userParams, bntParams)
    if ~isstruct(userParams)
        error('Incorrect argument');
    end

    bntSettingsFile = fullfile(helpers.bntRoot(), 'examples', 'bntSettings.m');
    if nargin < 2
        bntParams = io.readParamFile(bntSettingsFile);
    end

    bntFields = fieldnames(bntParams);
    userFields = fieldnames(userParams);

    p = userParams;

    missedFields = setdiff(bntFields, userFields);
    additionalFields = setdiff(userFields, bntFields);

    if ~isempty(missedFields)
        errMsg = sprintf('Your parameters file misses some values, defaults will be used. Here is the list of missing values:\n');
        errMsg = [errMsg sprintf('%s\n', missedFields{:})];
        warning('BNT:paramsDiff', [errMsg 'Check file %s for meaning of these fields.'], bntSettingsFile);
        for i = 1:length(missedFields)
            p.(missedFields{i}) = bntParams.(missedFields{i});
        end
    end
    if ~isempty(additionalFields);
        errMsg = sprintf('Your parameters file has some values, which are no longer in use. Here is the list of these values:\n');
        errMsg = [errMsg sprintf('%s\n', additionalFields{:})];
        warning('BNT:paramsDiff', [errMsg 'Check file %s for meaning of these fields.'], bntSettingsFile);
    end

    % let's check shuffle table which is complex and can be partially valid
    if ~isempty(find(strcmpi(userFields, 'shuffleTable'), 1))
        missedFields = setdiff(bntParams.shuffleTable(:, 1), userParams.shuffleTable(:, 1));
        additionalFields = setdiff(userParams.shuffleTable(:, 1), bntParams.shuffleTable(:, 1));
        if ~isempty(missedFields)
            warning('BNT:paramsDiff', 'Parameter p.shuffleTable in your parameters file misses some lines (definitions of shuffle variables). Check this parameter against the one in bntSettings.m');
        end
        if ~isempty(additionalFields)
            warning('BNT:paramsDiff', 'Parameter p.shuffleTable in your parameters file has more lines (definitions of shuffle variables) than are currently in use. You can safely remove them.\nCheck this parameter against the one in bntSettings.m');
        end
    end
end
