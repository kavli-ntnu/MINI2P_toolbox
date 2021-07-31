function pts = trianglePoints(dt, idx)
    trIndices = dt.ConnectivityList(idx, :);
    pts = zeros(length(trIndices), 2);
    for i = 1:length(trIndices)
        pts(i, :) = dt.Points(trIndices(i), :);
    end
end