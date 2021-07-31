% Check if the provided cache file is valid
%
% When the functions of BNT changes in such a way that previously created cache could be invalid,
% it is desirable to detect this situation. The function checks provided datetime and outputs TRUE if
% it is valid, or false if it is invalid. A valid datetime is older than the last valid datetime.
% The last valid datetime is hard coded in this function.
%
%  USAGE
%   isValid = helpers.isCacheValid(dateOfCreation)
%   dateOfCreation      Datetime that should be checked, serial datetime number.
%   isValid             TRUE or FALSE.
%
function isValid = isCacheValid(dateOfCreation)
    lastValidDate = datenum(2017, 05, 10, 12, 35, 00); % year, month, day, hour, minutes, secs

    diff_days = (lastValidDate) - (dateOfCreation);
    isValid = diff_days < 0;
end