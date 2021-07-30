function gscore = gridnessScoreShuffled(aCorr, cFieldRadius, gridnessRadius, radii)

    if cFieldRadius == 0
        gscore = nan;
        return;
    end

    halfSize = ceil(size(aCorr)/2);
    half_height = halfSize(1);
    half_width = halfSize(2);
    aCorrRad = min(halfSize);
    aCorrSize = size(aCorr);

    % contourc is efficient if aCorr is normalized
    maxValue = max(max(aCorr));
    if maxValue ~= 1
        aCorr = aCorr / maxValue;
    end

    % Meshgrid for expanding circle
    [rr, cc] = meshgrid(1:size(aCorr, 2), 1:size(aCorr, 1));

    % Define iteration radius step size for the gridness score
%     lastRadius = nanmin([aCorrRad maxRadius+3]);
    startInd = nanmax([cFieldRadius gridnessRadius-radii(1)]);
    endInd = nanmin([aCorrRad gridnessRadius+radii(2)]);
    radSteps = startInd:endInd;
    radSteps(1) = [];
    numSteps = length(radSteps);

    GNS = zeros(numSteps, 2);
    rotCorr = zeros(1, 5);
    rotAngles_deg = 30*(1:5);

    rotatedCorr = cell(1, length(rotAngles_deg));
    for i = 1:length(rotAngles_deg)
        rotatedCorr{i} = imrotate(aCorr, rotAngles_deg(i), 'bilinear', 'crop');
    end

    mainCircle = sqrt((cc - half_height).^2 + (rr - half_width).^2);
    innerCircle = mainCircle > cFieldRadius;

    % Define expanding ring of autocorrellogram and do x30 correlations
    for i = 1:numSteps
        ind = (innerCircle & (mainCircle < radSteps(i)));
        tempCorr = reshape(aCorr(ind), 1, [])';
        for j = 1:5
            rotatedCircle = reshape(rotatedCorr{j}(ind), 1, [])';
            rotCorr(j) =  corr(tempCorr, rotatedCircle);
        end
        GNS(i, 1) = min(rotCorr([2, 4])) - max(rotCorr([1, 3, 5]));
        GNS(i, 2) = radSteps(i);
    end

    % Find the biggest gridness score and radius
    gscore = nanmean(GNS(:, 1));
%     numGridnessRadii = 3;
%     numStep = numSteps - numGridnessRadii;
%     if numStep < 1
%         numStep = 1;
%     end
%
%     if numStep == 1
%         gscore = nanmean(GNS(:, 1));
%         % radius = nanmean(gridnessArray(:,2));
%     else
%         meanGridnessArray = zeros(numStep, 1);
%         for ii = 1:numStep
%             meanGridnessArray(ii) = nanmean(GNS(ii:ii + numGridnessRadii-1, 1));
%         end
%
%         [gscore, gInd] = max(meanGridnessArray);
%         % radius = gridnessArray(gInd+(p.numGridnessRadii-1)/2, 2);
%     end
end