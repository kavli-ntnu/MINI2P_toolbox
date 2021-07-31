% Group (split) spike timestamps into groups presented by vector edges
%
% This is just a wrapper around histc function with a more meaningful name for neuroscientists.
%
%  USAGE
%   [spkIndices] = general.groupSpikes(spikes, edges)
%   spikes          Spike timestamps, normally Nx1 matrix.
%   edges           Vector of edges that define group boundaries.
%                   See also general.groupByTime
%   spkIndices      Matrix Nx1 with logical indices of group distribution.
%                   spkIndices(3) == 1 means, that spike number 3 corresponds
%                   to group 1.
%
function spkIndices = groupSpikes(spikes, edges)
    if nargin < 2
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help general.groupSpikes">general.groupSpikes</a>'' for details).');
    end

    [~, spkIndices] = histc(spikes, edges);
end
