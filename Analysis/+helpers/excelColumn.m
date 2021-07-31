% Get column number in Excel format
%
% Use this function to convert a column number to Excel format. Excel columns
% are based on letters. For example, column number 2 in Excel will be 'B'.
%
%  USAGE
%   col = helpers.excelColumn(columnIndex)
%   columnIndex         Integer. Index of a column.
%   col                 String that represents column index in Excel format.
%
function col = excelColumn(columnIndex)
    columnIndex = columnIndex - 1;

    if columnIndex >= 0 && columnIndex < 26
        col = char('A' + columnIndex);
    elseif columnIndex > 25
        col = [helpers.excelColumn(floor(columnIndex / 26)) helpers.excelColumn(mod(columnIndex, 26) + 1)];
    else
        error('Invalid Excel column %u', columnIndex);
    end
end
