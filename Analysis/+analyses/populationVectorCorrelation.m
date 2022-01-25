% Calculates correlation between population vectors
%
% Population vector (stack) consists of rate maps. Each stack is a 3D
% matrix. The first two dimensions correspond to dimensions of a rate
% map for a cell. Third dimension is number of cells.
%
% This function returns different results based on input parameter 'full':
% *) 'vector' pvCorr is a vector in this case. stack1 and stack2 are processed
%     row- or column-wise depending on value of parameter 'orientation'.
%     For example, column processing means that each column i of each rate map
%     in a stack, is combined together in one long vector over all cells. As a
%     result there is long vector from stack1 and a long vector from stack2.
%     pvCorr(i) is correlation result of these two vectors (it is a single number).
%
% *) 'off' pvCorr is a 2D matrix in this case. Stacks are processed column-wise.
%    Each row i of pvCorr (pvCorr(i, :)) is a diagonal of a correlation matrix of
%    column i of stack1 and stack2.
%
% *) 'on' pvCorr is a 3D matrix in this case. Stacks are processed column-wise.
%    pvCorr consists of correlation matrices, which are created by correlation
%    columns of stack1 vs stack2. pvCorr(:, :, i) is a correlation matrix of
%    column i in stack1 and stack2.
%
%  USAGE
%   pvCorr = analyses.populationVectorCorrelation(stack1, stack2, <options>)
%   stack1      3D matrix, population vector
%   stack2      3D matrix, population vector
%   <options>   optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties    Values
%   -------------------------------------------------------------------------
%   'full'          'off' returns only diagonal elements of correlation matrix
%                   (i.e. correlation coefficients). pvCorr is a 2D matrix in
%                   this case. 'on' returns full correlation matrices. pvCorr
%                   is a 3D matrix in this case. 'vector' returns a vector:
%                   rate maps are processed row- or column-wise depending on
%                   value 'orientation'. Default is 'off'. (see description above)
%
%   'orientation'   'v' or 'h'. This parameter can only be used when 'full' is set to 'vector'.
%                   The parameter specifies how stacks are processed. 'v' means
%                   iteration over y-axis (rows) and calculation of correlation per row.
%                   'h' means iteration over x-axis (columns) and calculation of correlation
%                   per column. Default is 'v'.
%   'rows'          Rows parameter of Matlab's corr function. Refer to Matlab help for possible
%                   values. Default value is 'all'.
%   =========================================================================
%
%   pvCorr          Correlation values. Can be either 1D, 2D or 3D matrix based on
%                   value of parameter 'full'
%
function pvCorr = populationVectorCorrelation(stack1, stack2, varargin)

    if nargin < 2 || mod(length(varargin), 2) ~= 0,
      error('Incorrect number of parameters (type ''help <a href="matlab:help analyses.populationVectorCorrelation">analyses.populationVectorCorrelation</a>'' for details).');
    end

    returnFormats.D2 = 0;
    returnFormats.D3 = 1;
    returnFormats.D1 = 2;

    inp = inputParser;
    defaultFull = 'off';
    defaultOrientation = 'h';
    defaultRows = 'all';

    returnFormat = returnFormats.D2; % 0 - 2D matrix, 1 - 3D matrix, 2 - 1D vector
    isOrientationX = false;

    addRequired(inp, 'stack1');
    addRequired(inp, 'stack2');
    addParameter(inp, 'full', defaultFull, @(x) any(strcmpi({'on', 'off', 'vector'}, x)));
    addParameter(inp, 'orientation', defaultOrientation, @(x) any(strcmpi({'h', 'v'}, x)));
    addParameter(inp, 'rows', defaultRows, @ischar);

    parse(inp, stack1, stack2, varargin{:});
    if strcmpi(inp.Results.full, 'on')
        returnFormat = returnFormats.D3;
    elseif strcmpi(inp.Results.full, 'vector')
        returnFormat = returnFormats.D1;
    end
    if returnFormat == returnFormats.D1
        isOrientationX = strcmpi(inp.Results.orientation, 'h');
    end

    [yBins1, xBins1, numCells1] = size(stack1);
    [yBins2, xBins2, numCells2] = size(stack2);

    numCells = numCells1;
    if numCells1 ~= numCells2
        warning('BNT:stackSize', 'Population vectors are of different size: %u vs %u', numCells1, numCells2);
        numCells = min([numCells1 numCells2]);
    end

    numxBins = min([xBins1, xBins2]);
    numyBins = min([yBins1, yBins2]);

    switch returnFormat
        case returnFormats.D1
            if isOrientationX
                pvCorr = zeros(1, numxBins);
            else
                pvCorr = zeros(numyBins, 1);
            end

        case returnFormats.D2
            pvCorr = zeros(numyBins, numxBins);

        case returnFormats.D3
            pvCorr = zeros(xBins1, xBins2, numyBins);
    end

    if isOrientationX
        numBins = numxBins;
    else
        numBins = numyBins;
    end

    for i = 1:numBins % it doesn't matter whether we iterate over y or x bins
        if isOrientationX
            stack1dLeft = stack1(:, i, 1:numCells); % column
            stack1dRight = stack2(:, i, 1:numCells);
        else
            stack1dLeft = stack1(i, :, 1:numCells); % row
            stack1dRight = stack2(i, :, 1:numCells);
        end

        if returnFormat == returnFormats.D1
            if isOrientationX
                reshapedLeft = reshape(stack1dLeft, numCells * size(stack1dLeft, 1), 1)';
                reshapedRight = reshape(stack1dRight, numCells * size(stack1dRight, 1), 1)';
            else
                reshapedLeft = reshape(stack1dLeft, numCells * size(stack1dLeft, 2), 1)';
                reshapedRight = reshape(stack1dRight, numCells * size(stack1dRight, 2), 1)';
            end
        else
            % corr calculates correlation column-wise, so reshape
            % stacks and make them of size <numCells> x <number of bins>
            reshapedLeft = reshape(stack1dLeft, [], numCells)';
            reshapedRight = reshape(stack1dRight, [], numCells)';
        end

        switch returnFormat
            case returnFormats.D3
                pvCorr(:, :, i) = corr(reshapedLeft, reshapedRight, 'rows', inp.Results.rows);

            case returnFormats.D2
                pvCorr(i, :) = diag(corr(reshapedLeft, reshapedRight, 'rows', inp.Results.rows));

                % assign 0 to all NaN values
                pvCorr(i, isnan(pvCorr(i, :))) = 0;

            case returnFormats.D1
                if isOrientationX
                    pvCorr(1, i) = corr(reshapedLeft', reshapedRight', 'rows', inp.Results.rows);
                else
                    pvCorr(i, 1) = corr(reshapedLeft', reshapedRight', 'rows', inp.Results.rows);
                end
        end
    end
end
