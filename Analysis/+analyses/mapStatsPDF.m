% Calculate statistics of a rate map that depend on probability distribution function (PDF)
%
% Calculates information, sparsity and selectivity of a rate map. Calculations are done
% according to 1993 Skaggs et al. "An Information-Theoretic Approach to Deciphering the Hippocampal Code"
% paper. Another source of information is 1996 Skaggs et al. paper called
% "Theta phase precession in hippocampal neuronal populations and the compression of temporal sequences".
%
%  USAGE
%   [information, sparsity, selectivity] = analyses.mapStatsPDF(map)
%   map             Structure with rate map, output of analyses.map
%   information     Structure with fields:
%       rate        information rate [bits/sec]
%       content     Spatial information content [bits/spike]
%   sparsity        Sparsity value
%   selectivity     Selectivity value
%
%  SEE
%   See also analyses.map
%
function [information, sparsity, selectivity] = mapStatsPDF(map)
    if ~isstruct(map)
        error('BNT:arg', 'Incorrect argument. Map should be a structure. You are probably relying on old code when mapStatsPDF accepted map as matrix.');
    end
    T = nansum(nansum(map.time)); % overall trial duration
    posPDF = map.time / (T+eps); % probability of animal being in bin x:
                                 % duration of time animal spent in the bin divided by overall trial duration.
    
    meanrate = nansum(nansum(map.z .* posPDF));
    meansquarerate = nansum(nansum( (map.z .^ 2) .* posPDF ));
    if meansquarerate == 0
       sparsity = NaN;
    else
        sparsity = meanrate^2 / meansquarerate;
    end

    maxrate = nanmax(nanmax(map.z));
    if meanrate == 0;
       selectivity = NaN;
       information.content = nan;
       information.rate = nan;
    else
       selectivity = maxrate / meanrate;
       logArg = map.z / meanrate;
       logArg(logArg < 1) = 1;
       
       information.rate = nansum(nansum(posPDF .* map.z .* log2(logArg)));
       information.content = information.rate / meanrate;
    end
end
