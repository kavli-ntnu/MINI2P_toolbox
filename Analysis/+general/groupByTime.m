% Group (split) position data by time
%
% This function divides position data into groups based on value of groupParam parameter.
% You can use it to split position data into halves, or into groups of equal time.
%
%  USAGE
%   [groupIndices, edges] = general.groupByTime(pos, groupParam)
%   pos             Position samples. Matrix of size at least Nx3.
%   groupParam      Parameter that defines how groups are formed. Can be string
%                   or digit. The only possible value for string is 'half', which
%                   means that data will be divided into two groups.
%                   Numerical value defines group interval in seconds.
%   groupIndices    Vector of length equal to number of position samples. Each value
%                   indicates group belongingness. If groupIndices(3) == 1, means that
%                   position sample number 3 belongs to group number 1.
%   edges           Vector of edges that define groups. For example, if data is split 
%                   in halves, then edges = [0 <middle_point> Inf].
%
%  EXAMPLE
%
%   pos = data.getPositions();
%   groupIndices = general.groupByTime(pos, 'half'); % divide data into halves
%   positionsOfFirstGroup = pos(groupIndices == 1, :);
%
%   general.groupByTime(pos, 1); % divide data into groups each duration of 1 second
%
function [groupIndices, edges] = groupByTime(pos, groupParam)
    if nargin < 2
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help general.groupByTime">general.groupByTime</a>'' for details).');
    end

    if size(pos, 2) < 3
        error('BNT:arg', 'Incorrect argument ''pos'' (type ''help <a href="matlab:help general.groupByTime">general.groupByTime</a>'' for details).');
    end

    duration = pos(end, 1) - pos(1, 1);

    if ischar(groupParam)
        if ~strcmpi(groupParam, 'half')
            error('Incorrect argument ''groupParam'' (type ''help <a href="matlab:help general.groupByTime">general.groupByTime</a>'' for details).');
        end
        midPoint = pos(1, 1) + round(duration / 2);
        edges = [0 midPoint inf];
    else
        grpInterval = groupParam;
        if ~helpers.isdscalar(grpInterval, '>0')
            error('Incorrect argument ''groupParam'' (type ''help <a href="matlab:help general.groupByTime">general.groupByTime</a>'' for details).');
        end
        edges = [0:grpInterval:duration];
        edges(end+1) = inf;
    end

    [~, groupIndices] = histc(pos(:, 1), edges);
end
