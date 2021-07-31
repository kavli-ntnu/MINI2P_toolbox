% Identify 'bad' position samples based on animal's speed
%
% This scripts returns indices of position samples that are lower than lowSpeedThreshold
% or/and higher than highSpeedThreshold.
%
%  USAGE
%   ind = general.speedThreshold(pos, lowSpeedThreshold, highSpeedThreshold)
%   pos                     Position samples. Matrix of size at least Nx2.
%   lowSpeedThreshold       Lower bound of the speed threshold. Animal's speed should be greater 
%                           than this value. Provide 0 value if you do not want this part 
%                           of filtering.
%   highSpeedThreshold      Upper bound of the speed threshold. Animal's speed should be lower than
%                           this value. Provide 0 value if you do not want this part of filtering.
%   ind                     Indices of filtered samples, i.e. pos(ind, :) - are samples that should
%                           be filtered out.
%
function selected = speedThreshold(posIn, lowSpeedThreshold, highSpeedThreshold)
    if nargin < 3
        error('BNT:numArgs', 'Incorrect number of arguments');
    end
    if size(posIn, 2) < 2
        error('BNT:arg', 'Dimension of position data is incorrect. It should be a matrix of size Nx2.');
    end

    pos = posIn;
    if lowSpeedThreshold > 0 || highSpeedThreshold > 0
        v = general.speed(pos);
        if lowSpeedThreshold > 0 && highSpeedThreshold > 0
            selected = find(v < lowSpeedThreshold | v > highSpeedThreshold);
        elseif lowSpeedThreshold > 0
            selected = find(v < lowSpeedThreshold);
        else
            selected = find(v > highSpeedThreshold);
        end
    end
end
