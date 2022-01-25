%
% Returns indices of matched fields. These indices do not correspond to each other, i.e.
% there is no guarantee that matchedFields1(1) matches excatly matchedFields2(1). It matcehs
% one of matchedFields2.
%
function [matchedFields1, matchedFields2, matchedTriangles1, matchedTriangles2] = compareMapsByTriangulation(map1, fields1, map2, fields2)
    if isempty(fields1) || isempty(fields2)
        matchedFields1 = [];
        matchedFields2 = [];
        return;
    end
    
    dt1 = general.triangulateFields(fields1);
    dt2 = general.triangulateFields(fields2);

    numTriangles1 = size(dt1, 1);
    numTriangles2 = size(dt2, 1);
    distances = zeros(numTriangles2, 1);
    centroids2 = zeros(numTriangles2, 2);

    featDisHausdorf = zeros(numTriangles2, 1);
    featDistCityblock = zeros(numTriangles2, 1);

    matchedPoints1 = zeros(numTriangles1, 2);
    matchedPoints2 = zeros(numTriangles1, 2);
    matchedTriangles1 = cell(numTriangles1, 1); % coordinates of matched triangles in map1
    matchedTriangles2 = cell(numTriangles1, 1); % coordinates of triangles in map2 that match triangles in map1
    toRemove = false(numTriangles1, 1); % above matched data will be truncated based on toRemove variable

    matchedIndices = zeros(numTriangles1, 1);

    for i = 1:numTriangles1
        triangle1 = general.trianglePoints(dt1, i);
        matchedPoints1(i, :) = centroid(triangle1);
        matchedTriangles1{i} = triangle1;
        features1 = analyses.triangleFeatures(triangle1);

        %features2 = zeros(lengthfeatures1{i}, numTriangles2);

        for j = 1:numTriangles2
            triangle2 = general.trianglePoints(dt2, j);
            distances(j) = hausdorffDistance(triangle1, triangle2);
            centroids2(j, :) = centroid(triangle2);

            features2 = analyses.triangleFeatures(triangle2);
            featDisHausdorf(j) = hausdorffDistance(features1, features2);
            featDistCityblock(j) = pdist2(features1', features2', 'cityblock');
        end
        [minDist, minInd] = min(distances);
        if minDist < 6 && ~any(matchedIndices == minInd)
            matchedPoints2(i, :) = centroids2(minInd, :);
            matchedTriangles2{i} = general.trianglePoints(dt2, minInd);
            matchedIndices(i) = minInd;
        else
            toRemove(i) = true;
        end

%         [minDistHaus, minIndHaus] = min(featDisHausdorf);
%         fprintf('Iteration %u, min dist %f\n', i, minDistHaus);
%         if minDistHaus < 10
%             attempts = 1;
%             while any(matchedIndices == minIndHaus) && attempts < length(featDisHausdorf)
%                 featDisHausdorf(minIndHaus) = max(featDisHausdorf);
%                 attempts = attempts + 1;
%                 [minDistHaus, minIndHaus] = min(featDisHausdorf);
%             end
%             matchedPoints2(i, :) = centroids2(minIndHaus, :);
%             matchedTriangles2{i} = general.trianglePoints(dt2, minIndHaus);
%             matchedIndices(i) = minIndHaus;
%         else
%             toRemove(i) = true;
%         end
        
%         [minDistCity, minIndCity] = min(featDistCityblock);
%         if minDistCity > 80 && minDistCity < 100
%             attempts = 1;
%             while any(matchedIndices == minIndCity) && attempts < length(featDistCityblock)
%                 featDistCityblock(minIndCity) = max(featDistCityblock);
%                 attempts = attempts + 1;
%                 [minDistCity, minIndCity] = min(featDistCityblock);
%             end
%             matchedPoints2(i, :) = centroids2(minIndCity, :);
%             matchedTriangles2{i} = general.trianglePoints(dt2, minIndCity);
%             matchedIndices(i) = minIndCity;
%         else
%             toRemove(i) = true;
%         end
        
%         minIndCity = find(featDistCityblock > 480 & featDistCityblock < 500, 1);
%         if ~isempty(minIndCity)
%             attempts = 1;
%             while any(matchedIndices == minIndCity) && attempts < length(featDistCityblock)
%                 featDistCityblock(minIndCity) = max(featDistCityblock);
%                 attempts = attempts + 1;
%                 %[minDistCity, minIndCity] = min(featDistCityblock);
%                 minIndCity = find(featDistCityblock > 480 & featDistCityblock < 500, 1);
%                 if isempty(minIndCity)
%                     break;
%                 end
%             end
%             if ~isempty(minIndCity)
%                 matchedPoints2(i, :) = centroids2(minIndCity, :);
%                 matchedTriangles2{i} = general.trianglePoints(dt2, minIndCity);
%                 matchedIndices(i) = minIndCity;
%             else
%                 toRemove(i) = true;
%             end
%         else
%             toRemove(i) = true;
%         end
%         fprintf('Hausdorf %f, city %f, decision based on centroids: %u\n', minDistHaus, minDistCity, toRemove(i));
    end
    numLeft = numTriangles1 - length(find(toRemove));
%     fprintf('Found %u matched triangles, this is %f %%\n', numLeft, 
    matchedPoints1(toRemove, :) = [];
    matchedPoints2(toRemove, :) = [];
    matchedIndices(toRemove) = [];
    
    matchedTriangles1(toRemove) = [];
    matchedTriangles2(toRemove) = [];

    matchedFields1 = [];
    matchedFields2 = [];
    if ~isempty(matchedPoints1)
        matchedFields1 = matchedTriangles2Fields(find(~toRemove), dt1, fields1);
        matchedFields2 = matchedTriangles2Fields(matchedIndices, dt2, fields2);
%         figure;
%         plot.matchedTriangles(map1, map2, matchedPoints1, matchedPoints2, matchedTriangles1, matchedTriangles2);
    end
end

function fieldInd = matchedTriangles2Fields(indices, dt, fields)
    fieldInd = zeros(length(indices)*3, 1);
    fieldsX = [fields(:).x];
    fieldsY = [fields(:).y];
    
    numFields = 1;
    
    for i = 1:length(indices)
        triangle = general.trianglePoints(dt, indices(i));
        for j = 1:3
            candidates = (fieldsX == triangle(j, 1)) & (fieldsY == triangle(j, 2));
            ind = find(candidates, 1);
            if ~isempty(ind)
                fieldInd(numFields) = ind;
                numFields = numFields + 1;
            end
        end
    end
    
    fieldInd = unique(fieldInd, 'stable');
end