% Calculate 2D autocorrelation
%
% This script calculates 2D autocorrelation (autocorrelogram) of a firing map.
%
%  USAGE
%   Rxx = analyses.autocorrelation(map)
%   map         2D firing map
%   Rxx         Resulting correlation matrix
%
function Rxx = autocorrelation(map)
    overlapAmount = 0.8; % percentage of overlap region. We are interested only in
                         % overlapAmount of a full autocorrelogram. Edge areas produce artefacts.
                         % Shoulmapd be a value in range 0..1.

    map(isnan(map)) = 0;
    newSizeV = round(size(map, 1) + size(map, 1) * overlapAmount);
    newSizeH = round(size(map, 2) + size(map, 2) * overlapAmount);
    if mod(newSizeV, 2) == 0 && newSizeV > 0
        newSizeV = newSizeV - 1;
    end
    if mod(newSizeH, 2) == 0 && newSizeH > 0
        newSizeH = newSizeH - 1;
    end
    
    if isempty(map) || all(map(:) == 0)
        Rxx = zeros(newSizeV, newSizeH);
        return;
    end
    Rxx = SpatialTuning_BNT.normxcorr2_general(map, map);

    offsetV = size(Rxx, 1) -  newSizeV;
    offsetV = round(offsetV/2 + 1);
    offsetH = size(Rxx, 2) -  newSizeH;
    offsetH = round(offsetH/2 + 1);

    Rxx = Rxx(offsetV:end-offsetV+1, offsetH:end-offsetH+1); % +1 because offset should be taken
                                            % from size(map)*2, but it is taken from
                                            % size(map)*2 - 1
end
