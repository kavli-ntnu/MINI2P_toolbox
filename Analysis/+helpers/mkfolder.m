% Create folder if it doesn't exist
%
function mkfolder(folder)
    if exist(folder, 'dir') == 0
        mkdir(folder);
        if exist(folder, 'dir') == 0
            error('Failed to create output directory. Check the argument or create the directory manually');
        end
    end
end
