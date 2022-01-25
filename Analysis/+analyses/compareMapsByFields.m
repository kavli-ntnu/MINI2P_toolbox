function compareMapsByFields(map1, fields1, map2, fields2)
    % create features based on field properties, compare featrues
    
    numFields1 = length(fields1);
    numFields2 = length(fields2);
    distances = zeros(numFields2, 1);
    
    matchedPoints1 = zeros(numFields1, 2);
    matchedPoints2 = zeros(numFields1, 2);
    toRemove = false(numFields1, 1); % above matched data will be truncated based on toRemove variable
    
    matchedIndices = zeros(numFields1, 1);

    for i = 1:numFields1
        matchedPoints1(i, :) = [fields1(i).peakX fields1(i).peakY];
        features1 = fieldFeatures(fields1(i));
        
        for j = 1:numFields2
            features2 = fieldFeatures(fields2(j));
            distances(j) = hausdorffDistance(features1, features2);
%             peakPoint = [fields2(j).peakX fields2(j).peakY];
%             distances(j) = pdist2(matchedPoints1(i, :), peakPoint);
        end
        
        [minDist, minInd] = min(distances);
        if minDist < 21 && ~any(matchedIndices == minInd)
%         if ~any(matchedIndices == minInd)
            features2 = fieldFeatures(fields2(minInd));
            
            distance = pdist2(features1, features2);
            
            matchedPoints2(i, :) = [fields2(minInd).peakX fields2(minInd).peakY];
            matchedIndices(i) = minInd;
        else
            toRemove(i) = true;
        end
        fprintf('Field %u, distance %f, matched index %u\n', i, minDist, minInd);
        
    end
    matchedPoints1(toRemove, :) = [];
    matchedPoints2(toRemove, :) = [];

    if ~isempty(matchedPoints1)
        figure;
        plot.matchedPoints(map1, map2, matchedPoints1, matchedPoints2);
    else
        fprintf('No matches\n');
    end
end

function features = fieldFeatures(field)
%     features = zeros(5, 1);

    features(1) = field.x;
    features(2) = field.y;
%     features(3) = field.x;
%     features(4) = field.y;
%     features(3) = field.area;
%     features(1) = field.size;
    features(3) = field.Perimeter;
    features(4) = field.meanRate;
%     features(3) = field.Eccentricity;
%     features(4) = field.Extent;
%     features(5) = field.area;
    features(5) = field.Orientation;
end