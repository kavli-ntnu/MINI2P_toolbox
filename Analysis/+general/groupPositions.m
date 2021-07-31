% Group (split) position data by time or in specified number of partitions
%
% This function divides position data into groups based on pair of values (groupType, groupingValue).
% You can use it to split position data into predefined number of groups (2, 4, e.t.c.), or into
% groups of equal time.
%
%  USAGE
%   [groupIndices, edges] = general.groupPositions(pos, groupType, groupingValue)
%   pos             Position samples. Matrix of size at least Nx1.
%   groupType       String that defines type of grouping. Affects how groupingValue
%                   argument is used. 'time' defines group by time interval.
%                   'num' groups positions in specified number of groups.
%   groupingValue   Integer that specify value for grouping. If groupType == 'time'
%                   this argument should be a value of seconds.
%                   If groupType == 'num' this argument is the number of groups.
%
%   groupIndices    Vector of length equal to number of position samples. Each value
%                   indicates group belongingness. If groupIndices(3) == 1, means that
%                   position sample number 3 belongs to group number 1.
%   edges           Vector of edges that define groups. For example, if data is split 
%                   in halves, then edges = [0 <middle_point> Inf].
%
%  EXAMPLE
%
%   pos = data.getPositions();
%   groupIndices = general.groupPositions(pos, 'num', 2); % divide data into halves
%   positionsOfFirstGroup = pos(groupIndices == 1, :);
%
%   general.groupPositions(pos, 'time', 1); % divide data into groups each duration of 1 second
%
function [groupIndices, edges] = groupPositions(pos, groupType, groupingValue)
    if nargin < 3
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help general.groupPositions">general.groupPositions</a>'' for details).');
    end

    if size(pos, 1) == 1 && size(pos, 2) > 2
        % this is probably a row vector, transpose it
        pos = pos';
    end
    
    if size(pos, 2) < 1
        error('Incorrect value for argument ''pos'' (type ''help <a href="matlab:help general.groupPositions">general.groupPositions</a>'' for details).');
    end

    if ~helpers.isstring(groupType, 'time', 'num')
        error('Incorrect value for argument ''groupType'' (type ''help <a href="matlab:help general.groupPositions">general.groupPositions</a>'' for details).');
    end

    if ~helpers.isdscalar(groupingValue, '>0')
        error('Incorrect value for argument ''groupingValue'', it should be >= 1 (type ''help <a href="matlab:help general.groupPositions">general.groupPositions</a>'' for details).');
    end

    post = pos(:, 1);
    duration = max(post) - min(post);

    if strcmpi(groupType, 'num')
        interval = duration / groupingValue;
        edges = 0:interval:duration;
        edges = min(post) + edges; % shift if for case min(post) ~= 0
        edges(end) = Inf;
    else
        interval = groupingValue;
        edges = 0:interval:duration;
        edges = min(post) + edges; % shift if for case min(post) ~= 0
        edges(end+1) = Inf;
    end

    [~, groupIndices] = histc(post, edges);
end
