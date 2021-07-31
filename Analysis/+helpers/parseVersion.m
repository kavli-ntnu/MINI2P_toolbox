function [major, minor, numCommits, sha] = parseVersion(ver)
    if ~isempty(strfind(lower(ver), 'fatal'))
        error('BNT:arg', 'Incorrect argument ''ver'' (type ''help <a href="matlab:help helpers.parseVersion">helpers.parseVersion</a>'' for details).');
    end
    hyphenInd = strfind(ver, '-'); % should be two
    if isempty(hyphenInd)
        nums = regexpi(ver, '\d+', 'match'); % check that we have some numbers (should be two of them)
        if isempty(nums)
            error('BNT:arg', 'Incorrect argument ''ver'' (type ''help <a href="matlab:help helpers.parseVersion">helpers.parseVersion</a>'' for details).');
        end
        major = str2double(nums{1});
        minor = str2double(nums{2});
        numCommits = '';
        sha = '';        
    else
        if length(hyphenInd) ~= 2
            error('BNT:arg', 'Incorrect argument ''ver'' (type ''help <a href="matlab:help helpers.parseVersion">helpers.parseVersion</a>'' for details).');
        end
        verId = ver(1:hyphenInd(1)-1);
        numCommits = ver(hyphenInd(1)+1:hyphenInd(2)-1);
        sha = ver(hyphenInd(2)+2:end); % remove 'g' character
        nums = regexp(verId, '\d*', 'match');
        major = str2double(nums{1});
        minor = str2double(nums{2});
    end
end