% Identify trial by cell number
%
% This function searches for a trial to which a cell belongs. Cell is
% identified by it's linear index. The function prints out information
% about the trial.
%
%  USAGE
%   helpers.cell2trial(cell)
%   cell    Cell number, unsigned integer.
%
function cell2trial(cell)
    if cell <= 0
        fprintf('Cell %i was not found within loaded trials\n', cell);
        return;
    end
    N = 0;
    target = cell;
    for i = 1:data.numTrials
        data.setTrial(i);
        units = data.getCells();
        if target > N && target <= N+size(units, 1)
            fprintf('Cell %u belongs to trial %u\n', target, i);
            data.getCurrentSession
            return;
        end
        N = N + size(units, 1);    
    end
    fprintf('Cell %u was not found within loaded trials\n', target);
end