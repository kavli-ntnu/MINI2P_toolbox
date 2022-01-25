% Calculate adaptive smoothed rate map
%
% Calculates an adaptive smoothed rate map as described in "Skaggs et al 1996 -
% Theta Phase Precession in Hippocampal Neuronal Population and the Compression of Temporal Sequences"
%
%  USAGE
%   [map, posPdf] = analyses.ratemapAdaptiveSmoothing(pos, spkPos, <options>)
%   pos         Position samples (either Nx3 or Nx2 matrix). X and Y coordinates are used.
%   spkPos      Spike position samples (either Nx2 or Nx2 matrix). X and Y coordinates are used.
%   <options>   optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties     Values
%   -------------------------------------------------------------------------
%    'binWidth'     width (and height) of bins in firing rate map. Value units
%                   are [cm]. binWidth has priority over parameter nBins. See
%                   description of nBins.
%
%    'nBins'        number of horizontal and vertical bins (default = [50 50]).
%                   If single value is provided, it is used for both horizontal
%                   and vertical number of bins (i.e. nBins = 50 => [50 50]).
%                   Either binWidth or nBins should be provided. If non is provided,
%                   then the default value of binWidth is used. If both are provided,
%                   then the value of binWidth is used.
%
%    'minTime'      minimum time spent in each bin (in s, default = 0).
%    'alphaValue'   scaling parameter. Default value is 10000. See original paper
%                   for the details.
%    'shape'        shape of the arena. Possible values are: 1 for square box,
%                   2 for cylinder. Only square box is currently supported!
%    =========================================================================
%
%  OUTPUT
%
%   map.x       x bins
%   map.y       y bins
%   map.z       adaptively smoothed firing rate map
%   posPdf      position probability density function

function [map, posPdf] = mapAdaptiveSmoothing(pos, spkPos, varargin)
    % Check number of parameters
    if nargin < 2 || mod(length(varargin), 2) ~= 0,
        error('BNT:numArgs', 'Incorrect number of parameters (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
    end

    % Check parameter sizes
    if size(pos, 2) < 2
        error('Parameter ''pos'' should have at least 2 columns (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
    end
    if size(spkPos, 2) < 2
        error('Parameter ''spkPos'' should have 1 or 2 columns (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
    end
    
    if isempty(pos)
        return;
    end

    % Default values
    alphaValue = 10000;
    binWidth = 5;
    minTime = 0;
    limits = [];
    shape = 1;
    sampleTime = helpers.sampleTimeFromData(pos);
    nBins = [50 50];

    haveBinWidth = false;
    haveNBins = false;

    % Parse parameter list
    i = 1;
    while i < length(varargin)
        if ~ischar(varargin{i}),
            error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).']);
        end

        switch(lower(varargin{i})),
            case 'binwidth',
                binWidth = varargin{i+1};
                if ~helpers.isdscalar(binWidth, '>0')
                    error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
                end
                haveBinWidth = true;
                i = i + 2;

            case 'nbins'
                nBins = varargin{i+1};
                if ~helpers.isivector(nBins, '>0') || length(nBins) > 2,
                    error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
                end
                haveNBins = true;
                i = i + 2;

            case 'mintime',
                minTime = varargin{i+1};
                if ~helpers.isdscalar(minTime, '>=0'),
                    error('Incorrect value for property ''minTime'' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
                end
                i = i + 2;

            case 'alphavalue'
                alphaValue = varargin{i+1};
                if ~helpers.isiscalar(alphaValue, '>0'),
                    error('Incorrect value for property ''alphaValue'' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
                end
                i = i + 2;
                
            case 'limits'
                limits = varargin{i+1};
                if ~isvector(limits)
                    error('Incorrect value for property ''limits'' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).');
                end
                i = i + 2;

            otherwise,
                error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help analyses.mapAdaptiveSmoothing">analyses.mapAdaptiveSmoothing</a>'' for details).']);
        end
    end

    % Some info about x, y and z
    if size(pos, 2) > 2
        posx = pos(:, 2);
        posy = pos(:, 3);
    else
        posx = pos(:, 1);
        posy = pos(:, 2);
    end

    if size(spkPos, 2) > 2
        spkx = spkPos(:, 2);
        spky = spkPos(:, 3);
    else
        spkx = spkPos(:, 1);
        spky = spkPos(:, 2);
    end

    if length(nBins) == 1
        nBins(2) = nBins(1);
    end

    % Calculate the border coordinates
    % xStart        Minimum x-coordinate for the path
    %
    % yStart        Minimum y-coordinate for the path
    %
    % xLength       Length of the arena in the x-direction [cm](for cylinder 
    %               this equals the diameter)
    % yLength       Length of the arena in the y-direction [cm] (for cylinder
    %               this equals the diameter)
    maxX = nanmax(posx);
    maxY = nanmax(posy);
    xStart = nanmin(posx);
    yStart = nanmin(posy);
    xLength = maxX - xStart;
    yLength = maxY - yStart;

    if haveNBins && ~haveBinWidth
        binWidth = ceil(xLength / nBins(1));
        numColBins = nBins(1);
        numRowBins = nBins(2);
    else
        % Number of bins in each direction of the map
        numColBins = ceil(xLength / binWidth);
        numRowBins = ceil(yLength / binWidth);
    end

    halfBin = binWidth / 2;

    rowAxis = halfBin:binWidth:(numRowBins*binWidth)-halfBin;
    rowAxis = yStart + rowAxis';
    colAxis = halfBin:binWidth:(numColBins*binWidth)-halfBin;
    colAxis = xStart + colAxis';

    maxBins = max([numColBins, numRowBins]);

    map.x = colAxis;
    map.y = rowAxis;
    map.z = zeros(numRowBins, numColBins);
    posPdf = zeros(numRowBins, numColBins);

    if shape(1) == 1
        % Overall clue:
        %     - grow circle from r=1:maxBins (mult. of binWidth), tracking inside
        %     - stop at smallest rad. such that r >= alpha/samples*sqrt(spikes)
        radsqs = ((1:maxBins) * binWidth) .^ 2; % square radius in cm
        
        for i = 1:numColBins
            binPosX = colAxis(i);
            
            dist_sample_xdir = (posx - binPosX) .^ 2;
            dist_spike_xdir = (spkx - binPosX) .^ 2;
            
            for j = 1:numRowBins
                binPosY = rowAxis(j);
                
                % Calculate sample and spike distances from bin center
                dist_sample = dist_sample_xdir + (posy-binPosY).^2;
                dist_spike = dist_spike_xdir + (spky-binPosY).^2;

                found = 0;
                % Grow circle in increments of binWidth
                for r = 1:maxBins
                   n = nnz(dist_sample <= radsqs(r));
                   s = nnz(dist_spike <= radsqs(r));
                   
                   if r >= alphaValue/(n*sqrt(s))
                       found = 1;
                       break;
                   end
                end
     
                % Set the rate for this bin
                map.z(j, i) = found * s/(n*sampleTime);
                posPdf(j, i) = found * n*sampleTime;
            end
        end 
    else
        for ii = 1:numColBins
            binPosY = (yStart + binWidth/2);
            for jj = 1:numRowBins
                currentPosition = sqrt(binPosX^2 + binPosY^2);
                if currentPosition > shape(2)/2
                    map.z(numRowBins-jj+1, ii) = NaN;
                    posPdf(numRowBins-jj+1, ii) = NaN;
                else
                    n = 0;
                    s = 0;
                    for r = 1:maxBins
                        % Set the current radius of the circle
                        radius = r * binWidth;
                        % Number of samples inside the circle
                        n = insideCircle(binPosX, binPosY, radius, posx, posy);
                        % Number of spikes inside the circle
                        s = insideCircle(binPosX, binPosY, radius, spkx, spky);

                        if r >= alphaValue/(n*sqrt(s))         
                            break;
                        end

                    end
                    % Set the rate for this bin
                    map.z(jj,ii) = s/(n*sampleTime);
                    posPdf(jj,ii) = n*sampleTime;
                    
                end
                binPosY = binPosY + binWidth;
            end 

            binPosX = binPosX + binWidth;
        end
    end

    map.z(posPdf < minTime) = NaN;
    posPdf = posPdf / nansum(nansum(posPdf));
end
