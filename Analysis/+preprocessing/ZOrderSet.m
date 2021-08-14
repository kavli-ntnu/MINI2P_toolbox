%ZOrderSet   Set a new z-order position of an existing object of a figure
%
%   ZOrderSet(h = gco, zIndex = 0 (last, default) | n (n-th) | -n (n-th last))
%
%   Set a new z-order position of an existing object of a figure.
%   The lower the z-order, the closer the object to the user eye.
%
%   The zIndex is the target z-order of the object. If out of boundaries,
%   it will be clipped.
%
%   When adding a new object to a figure, by default it is provided with
%   the first z-order value, so it becomes the most visible one.
%
%   Example:
%      figure(1); clf; hold on;
%      h1 = plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0]);
%      h2 = plot(1, 1, 'o', 'MarkerSize', 20, 'MarkerFaceColor', [0 1 0]);
%      % the red marker (h1) is not visible at this time
%      ZOrderSet(h1, 1); % set it to first position: the red marker is visible
%      ZOrderSet(h1, 2); % set it to 2nd (last) position: the red marker is not visible
%      ZOrderSet(h1, 0); % set it to last position: the red marker is not visible
%
%   Thanks and inspired to AddReorderButtons by Geoffrey K. Adams, November 2007
%
%   Copyright 2011
%
%   v1.0.0 - 28/04/2011
%   Marcello Ferro <marcello.ferro@ilc.cnr.it>
%   http://www.ilc.cnr.it
%
function ZOrderSet(h, zIndex)
% Check params
if(nargin < 2)
    zIndex = 0;
end
if(nargin < 1)
    h = gco;
end
% Get the parent of the object
parent = get(h, 'Parent');
if(parent == 0)
    % This means the object is probably the figure... So just ignore
    return;
end
% Get all children of the object's parent
children = get(parent, 'Children');
% Count children
hCount = length(children);
% Adjust target if negative
if(zIndex < 1)
    zIndex = hCount + zIndex;
end
% Clip target to boundaries
zIndex = min(zIndex, hCount);
zIndex = max(zIndex, 1);
% Find out this object's index number
hIndex = find(h == children);  
% If it is already first, it can't be moved forward any more
if(hIndex == zIndex)                       
    return;                        
end
if(zIndex < hIndex)
    % Move a slice forward to make place for the target
    children((zIndex+1):hIndex) = children(zIndex:(hIndex-1));
else
    % Move a slice backward to make place for the target
    children(hIndex:(zIndex-1)) = children((hIndex+1):zIndex);
end
% Place the target
children(zIndex) = h; 
% Set the new children
set(parent, 'Children', children);
end