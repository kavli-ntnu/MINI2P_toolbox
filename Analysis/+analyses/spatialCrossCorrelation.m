% Calculate Pearson cross correlation between two firing rate maps
%
% The size of return variable corrValue depends on input parameter 'output':
% *) 'single', corrValue is a single number. Maps are reshaped in a vector and
%     correlation is calculated.
% *) 'vector', corrValue is a vector. Maps are processed row- or column-wise.
%    corrValue(i) is a correlation coeeficient for row/column i.
% Row/column selection is done by parameter 'processBy'.
%
%  USAGE
%   corrValue = analyses.spatialCrossCorrelation(map1, map2, <options>)
%   map1      2D matrix, rate map
%   map2      2D matrix, rate map
%   <options>   optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties    Values
%   -------------------------------------------------------------------------
%   'output'       'single' or 'vector'. Default is 'single'.
%
%   'processBy'     'r' or 'c'. The parameter specifies how maps are processed.
%                   'r' means iteration over y-axis (rows) and calculation of correlation per row.
%                   'c' means iteration over x-axis (columns) and calculation of correlation
%                   per column. Default is 'c'. This only makes difference for 'vector' output.
%                   Single value is independent of processing order.
%   =========================================================================
%
%   corrValue       Correlation value(s). Can be either a single number or a vector.
%
function corrValue = spatialCrossCorrelation(map1, map2, varargin)
    columnWise = 0;
    rowWise = 1;
    single = 'single';

    inp = inputParser;
    defaultProcessBy = 'c';
    defaultReturnFormat = single;

    checkProcessBy = @(x) strcmpi(x, 'c') || strcmpi(x, 'r');
    checkReturnFormat = @(x) strcmpi(x, 'single') || strcmpi(x, 'vector');

    addRequired(inp, 'map1');
    addRequired(inp, 'map2');
    addParamValue(inp, 'processBy', defaultProcessBy, checkProcessBy);
    addParamValue(inp, 'output', defaultReturnFormat, checkReturnFormat);

    parse(inp, map1, map2, varargin{:});

    % get parsed arguments
    if strcmpi(inp.Results.processBy, 'c')
        orientation = columnWise;
    else
        orientation = rowWise;
    end
    returnFormat = inp.Results.output;

    if ~isequal(size(map1), size(map2))
        warning('BNT:mapsSize', 'Maps are of different size. This is suspicious. Correlation won''t be calculated');
        corrValue = nan;
        return;
    end

    if strcmpi(returnFormat, single)
        % Transform the 2-D maps to 1-D arrays by assembling the columns/rows from the maps
        A = reshape(map1, numel(map1), 1);
        B = reshape(map2, numel(map2), 1);

        % Remove the bins with NaN in A from both maps
        nansA = isnan(A);
        A = A(~nansA);
        B = B(~nansA);

        % Remove the bins with NaN in B from both maps
        nansB = isnan(B);
        A = A(~nansB);
        B = B(~nansB);

        % Calculate the correlation for the two rate maps.
        if isempty(A) || isempty(B)
            corrValue = NaN; % Cannot do calculation on empty arrays
        else
            corrValue = corr(A, B);
        end
    else
        nansLeft = isnan(map1);
        map1(nansLeft) = [];
        map2(nansLeft) = [];

        nansRight = isnan(map2);
        map1(nansRight) = [];
        map2(nansRight) = [];

        if orientation == columnWise
            corrMatrix = corr(map1, map2);
            corrValue = diag(corrMatrix)'; % make a row
        else
            corrMatrix = corr(map1', map2');
            corrValue = diag(corrMatrix); % make a column
        end
    end
end
