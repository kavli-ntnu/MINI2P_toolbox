function runs = findLaps(position)
    % create an empty array of tracks
    tracks = struct(...
        'id', {}, ...
        'kalmanFilter', {}, ...
        'age', {}, ...
        'totalVisibleCount', {}, ...
        'consecutiveInvisibleCount', {});
    posx = position(:, 2);
    posy = position(:, 3);
    numSamples = length(posx);
    
    fprintf('Finding individual runs. This takes some time...');

    minVisibleCount = 10;
    nextId = 1; % ID of the next track
    
    visualize = false;
    
    if visualize
        h = figure();
        axis([-110 110 -20 20]);
        hold on;
    end
    
    lastId = 0;
    endTrack = true;
    
    runs = zeros(0, 2);
    runInd = 1;
    
    for s = 1:numSamples
        % detect object
        centroids = [posx(s) posy(s)]; % 1x2 matrix
        if sum(isnan(centroids)) > 0
            continue;
        end
%         if visualize, plot(posx(s), posy(s), 'g*'), end

        predictNewLocations();
        [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment();

        updateAssignedTracks();
        updateUnassignedTracks();
        deleteLostTracks();
        createNewTracks();
        
        if ~isempty(tracks)
            reliableTrackInds = [tracks(:).totalVisibleCount] > minVisibleCount;
            reliableTracks = tracks(reliableTrackInds);

            if ~isempty(reliableTracks)
                ids = int32([reliableTracks(:).id]);
                if length(reliableTracks) > 1
                    fprintf('Several tracks!\n');
                end
                ids = ids(1);
                if lastId ~= ids
                    lastId = ids;
                    endTrack = true;
                    runs(runInd, 1) = s;
                    
%                     fprintf('Got new position index %d\n', s);
                    if visualize
                        cla
                    end
                end
                if visualize
                    plot(posx(s), posy(s), 'g*');
                    pause(0.001);
                end
            end
        end
        
%         if visualize, pause(0.001), end
    end
    
    if visualize, close(h), end
    
    if runs(runInd, 2) == 0
        runs(runInd, 2) = length(posx);
    end
    
    fprintf('done\n');

    function predictNewLocations()
        for i = 1:length(tracks)
            predictedCentroid = predict(tracks(i).kalmanFilter);
%             fprintf('Predicted %f\n', predictedCentroid);
%             if visualize
%                 plot(predictedCentroid(1), predictedCentroid(2), 'r+');
%             end
        end
    end

    function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment()
        
        nTracks = length(tracks);
        nDetections = size(centroids, 1);
        
        % compute the cost of assigning each detection to each track
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
        end
        
        % solve the assignment problem
        costOfNonAssignment = 10;
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);
    end

    function updateAssignedTracks()
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            centroid = centroids(detectionIdx, :);
            % bbox = bboxes(detectionIdx, :);
            
            % correct the estimate of the object's location
            % using the new detection
            correct(tracks(trackIdx).kalmanFilter, centroid);
            
            % update track's age
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            
            % update visibility
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            tracks(trackIdx).consecutiveInvisibleCount = 0;
        end
    end

    function updateUnassignedTracks()
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
        end
    end

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 10;
        ageThreshold = 8;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.6) | ...
            [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        if sum(lostInds) > 0 && endTrack
%             fprintf('End of track %d\n', s);
            if size(runs, 1) >= runInd
%                 if s - 8 < runs(runInd, 1)
%                     a = 5;
%                 end
                runs(runInd, 2) = s - 8;
                runInd = runInd + 1;
                endTrack = false;
            end
        end
        
        % delete lost tracks
        tracks = tracks(~lostInds);
    end

    function createNewTracks()
        centroids = centroids(unassignedDetections, :);
        
        for i = 1:size(centroids, 1)
            centroid = centroids(i,:);
            
            % create a Kalman filter object
            kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                centroid, [50, 15], ... % location variance, velocity variance
                [10, 5], 50);
            
            % create a new track
            newTrack = struct(...
                'id', nextId, ...
                'kalmanFilter', kalmanFilter, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0);
            
            % add it to the array of tracks
            tracks(end + 1) = newTrack;
            
            % increment the next id
            nextId = nextId + 1;
        end
    end
end
