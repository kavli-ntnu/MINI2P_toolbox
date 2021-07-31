% Find the centre of the box
function centre = centreBox(posx, posy)
    % Find border values for path and box
    maxX = max(posx);
    minX = min(posx);
    maxY = max(posy);
    minY = min(posy);
    
    centre = [nan nan];
    if maxX == minX || all(isnan([maxX minX]))
        centre(1) = 0;
    end
    if maxY == minY || all(isnan([maxY minY]))
        centre(2) = 0;
    end
    
    if centre(1) == 0 && centre(2) == 0
        return;
    elseif centre(1) == 0
        centre(2) = minY + (maxY - minY) / 2;
        return;
    elseif centre(2) == 0
        centre(1) = minX + (maxX - minX) / 2;
        return;
    end

    % Set the corners of the reference box
    NE = [maxX, maxY];
    NW = [minX, maxY];
    SW = [minX, minY];
    SE = [maxX, minY];

    % Get the centre coordinates of the box
    centre = findCentre(NE,NW,SW,SE);
end

% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE)

    % The centre will be at the point of interception by the corner diagonals
    a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
    b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
    c = SW(2);
    d = NW(2);
    x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
    y = a*(x-SW(1))+c; % Y-coord of centre
    centre = [x,y];
end
