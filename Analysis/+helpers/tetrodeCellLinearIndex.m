% Map tetrode and cell number to an integer
%
% Function produces a linear index from tetrode and cell number.
%
function ind = tetrodeCellLinearIndex(tetrode, cell)
    bigNumber = 1e6;
    ind = sub2ind([bigNumber bigNumber], tetrode, cell);
end