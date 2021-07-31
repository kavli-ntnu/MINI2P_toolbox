% Run an update procedure on a Git repository
%
% This function runs a Git update (pull) process on a given local repository.
% It was intended to be used in BNT, but can be used with other repositories.
% The update is carried out based on a `master` branch. This is hard coded.
%
%  VERSIONS
%   The whole update process is based on Git tag names that should contain version information.
%   A kind of verMajor.verMinor-patch schema is used for versioning, which resembles semantic versioning.
%   The patch version is not applied manually, but it is a number of commits since last version tag.
%   A valid tag for example is `v1.2`. From such a tag, major and minor versions are extracted.
%   Major is 1 and minor is 2 in this example. Minor and patch updates are carried out automatically.
%   For a major update, a user is asked for confirmation.
%
%  UPDATE INTERVALS AND FORCED UPDATES
%   By default, this function doesn't rerun an update if the last one has been carried out
%   less than 3 hours ago (time interval of 3 hours is hard coded). The function does this
%   by writing the time of last update in a file `repoName-updCheckTime.txt` in system temp
%   folder. The function then checks time from that file against current time and decides
%   whether update should continue.
%
%  NO UPDATES
%   It is possible to disable updates. The function check existence of file `no_update`
%   under repoPath. If there is one, then no update is performed. For BNT related projects
%   this file is checked under BNT root folder.
%
%  USAGE
%   helpers.runGitUpdate(repoName, repoPath, force)
%   repoName    Repository name, like 'BNT'. If not provided,
%               'BNT' is used.
%   repoPath    Full path to the repository or .git folder on hard drive.
%               If not provided, BNT root folder is used.
%   force       Flag that indicates whether we should force an update
%               or not. True is to force it. Default value is false.
%               See above about update intervals.
%   updated     Flag that indicated whether an update has been successful or not.
%
function updated = runGitUpdate(repoName, repoPath, force)
    if nargin < 1
        repoName = 'BNT';
    end

    if nargin < 2
        repoPath = helpers.bntRoot();
    end

    branchName = 'master';

    if nargin < 3
        force = false;
    end

    if ~ischar(repoName)
        error('Incorect parameter ''repoName''. It should be a string');
    end
    if ~ischar(branchName)
        error('Incorect parameter ''branchName''. It should be a string');
    end
    if ~ischar(repoPath)
        error('Incorect parameter ''repoPath''. It should be a string');
    end

    if isempty(strfind(repoPath, '.git'))
        repoGitPath = fullfile(repoPath, '.git');
    end

    updated = false;
    noUpdateInterval_sec = 10800; % 3 hours

    fprintf('Checking for updates: %s...\n', repoName);

    if ~isempty(strfind(lower(repoName), 'bnt'))
        noUpdatesFile = fullfile(helpers.bntRoot, 'no_updates');
    else
        noUpdatesFile = fullfile(repoPath, 'no_updates');
    end
    if exist(noUpdatesFile, 'file') ~= 0
        fprintf('Found ''no_updates'' file. Autoupdate will be stopped.\n');
        return;
    end

    % check when the last update took place
    timeFile = fullfile(tempdir, [repoName '-updCheckTime.txt']);
    if exist(timeFile, 'file') ~= 0
        fid = data.safefopen(timeFile, 'r');
        lastUpdateTime = fscanf(fid, '%f');
        clear fid;
        diffWithNow_sec = etime(datevec(now), datevec(lastUpdateTime));
        if force == false && diffWithNow_sec < noUpdateInterval_sec
            fprintf('Will skip the check. The last one have been performed recently\n');
            return;
        end
    end

    sl = git(sprintf('--git-dir "%s" --work-tree="%s" fetch', repoGitPath, repoPath)); %#ok<*NASGU>
    if exist(repoPath, 'dir') == 0 || ~isempty(strfind(sl, 'fatal'))
        % there is no .git folder, so we can not really proceed
        warning('BNT:updateFailed', 'Code seems to be not under version control. Can not run an automatic update in this case.');
        return;
    end

    branchName = git(sprintf('--git-dir "%s" --work-tree="%s" rev-parse --abbrev-ref HEAD', repoGitPath, repoPath));

    % get version information
    curVersion = git(sprintf('--git-dir "%s" --work-tree="%s" describe --long --tags --match v*', repoGitPath, repoPath));
    remoteVersion = git(sprintf('--git-dir "%s" --work-tree="%s" describe --long --tags --match v* origin/%s', repoGitPath, repoPath, branchName));
    masterVersion = git(sprintf('--git-dir "%s" --work-tree="%s" describe --long --tags --match v* origin/master', repoGitPath, repoPath));

    [curMajor, curMinor, curPassedCommits, curSha] = helpers.parseVersion(curVersion);
    [remMajor, remMinor, remPassedCommits, remSha] = helpers.parseVersion(remoteVersion); %#ok<*ASGLU>
    [masterMajor, masterMinor, masterPasstCommits, masterSha] = helpers.parseVersion(masterVersion);

    latestCommitSha = '';

    if curMajor ~= remMajor && strcmpi(branchName, 'master')
        % no automatic major updates on master. However if not on master a minor update can still be done

        % check for presence of minor updates before the major
        nextMajor = sprintf('v%u.0', curMajor + 1);
        latestCommit = git(sprintf('--git-dir "%s" --work-tree="%s" log --pretty=oneline --reverse -1 %s..%s~1', repoGitPath, repoPath, curSha, nextMajor));
        if ~isempty(latestCommit)
            latestCommitSha = strtok(latestCommit, ' ');
            % reinitialize remote version string to point to the latest minor update
            remoteVersion = git(sprintf('--git-dir "%s" --work-tree="%s" describe --long --tags --match v* %s', repoGitPath, repoPath, latestCommitSha));
            [remMajor, remMinor, remPassedCommits, remSha] = helpers.parseVersion(remoteVersion);
        end
    end
    if curMajor ~= masterMajor
        fprintf('There is a major update available. Do you want to update?\n');
        answ = input('Type ''y'' to update or any other letter NOT to do the update: ', 's');
        if strcmpi(answ, 'y')
            latestCommitSha = ''; % reset it, so we get the latest commit
        else
            if isempty(latestCommitSha)
                % only exit if there is no latestCommitSha. Otherwise do a minor update
                saveLastUpdateTime();
                return;
            end
        end
    end

    curPassedCommits = str2double(curPassedCommits);
    remPassedCommits = str2double(remPassedCommits);
    if curMinor == remMinor && curPassedCommits == remPassedCommits
        fprintf('Your version is the current one on selected branch!\n');
        saveLastUpdateTime();
        return;
    end

    doUpdate = true;

    diffStatus = git(sprintf('--git-dir "%s" --work-tree="%s" diff', repoGitPath, repoPath));
    if ~isempty(diffStatus)
        fprintf('You have some local changes. Update will erase these changes.\nDo you want to continue with the update?\n');
        answ = input('Type ''y'' to continue or any other letter to stop the update: ', 's');
        if ~strcmpi(answ, 'y')
            doUpdate = false;
        end
    end

    if ~doUpdate
        saveLastUpdateTime();
        return;
    end

    % clear our mex files
    clear Mat2Nlx* Nlx2Mat*;

    % revert all local changes
    sl = git(sprintf('--git-dir "%s" --work-tree="%s" checkout --force %s', repoGitPath, repoPath, branchName));

    sl = git(sprintf('--git-dir "%s" --work-tree="%s" pull --no-rebase', repoGitPath, repoPath)); % pull can easily jump to a next major version
    if ~isempty(latestCommitSha)
        sl1 = git(sprintf('--git-dir "%s" --work-tree="%s" reset --hard %s', repoGitPath, repoPath, latestCommitSha));
    end

    newVersion = git(sprintf('--git-dir "%s" --work-tree="%s" describe --long --tags --match v*', repoGitPath, repoPath));
    fprintf('Updated successfully from version %s to %s.\nIf you were running any script, rerun it!\n', curVersion, newVersion);
    updated = true;

    saveLastUpdateTime();

    %% Write time of last check in file
    function saveLastUpdateTime()
        try
            fid = data.safefopen(timeFile, 'w');
            fprintf(fid, '%f', now);
            clear fid;
        catch err
        end
    end
end
