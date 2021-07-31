% Fit square (rectangle) to provided points
function sq = fitSquare(points)
    goodInd = ~isnan(points(:, 1));
    points = points(goodInd, :);

    %% cluster points and arrange them counter-clockwise from closest to [0, 0]
    % rectangle:
    %   d    c
    %   a    b
    [clusterIdx, c] = kmeans(points, 4);
    d = pdist2(c, [0 0]);

    [~, idx] = min(d); % minimum distance from the origin
    aIdx = idx;

    selection = true(size(c, 1), 1);
    selection(idx) = false;

    % get lower-right corner
    [~, sIdx] = sort(c(:, 1), 'descend');
    [~, pointerToPointer] = min(c(sIdx(1:2), 2));
    idx = sIdx(pointerToPointer);
    selection(idx) = false;
    bIdx = idx;

    % get upper-right corner
    [~, sIdx] = sort(c(:, 2), 'descend');
    [~, pointerToPointer] = max(c(sIdx(1:2), 1));
    idx = sIdx(pointerToPointer);
    selection(idx) = false;
    cIdx = idx;

    dIdx = find(selection);

    %% now assemble point that lie on the same sides of rectangle
    % P = da, Q = ab, R = bc, S = cd
    idx = clusterIdx == dIdx | clusterIdx == aIdx;
    Px = points(idx, 1);
    Py = points(idx, 2);

    idx = clusterIdx == aIdx | clusterIdx == bIdx;
    Qx = points(idx, 1);
    Qy = points(idx, 2);

    idx = clusterIdx == bIdx | clusterIdx == cIdx;
    Rx = points(idx, 1);
    Ry = points(idx, 2);

    idx = clusterIdx == cIdx | clusterIdx == dIdx;
    Sx = points(idx, 1);
    Sy = points(idx, 2);

    %% fit rectangle
    zp = zeros(size(Px)); op =  ones(size(Px));
    zq = zeros(size(Qx)); oq =  ones(size(Qx));
    zr = zeros(size(Rx)); or =  ones(size(Rx));
    zs = zeros(size(Sx)); os =  ones(size(Sx));

    A = [ op zp zp zp Px  Py
             zq oq zq zq Qy -Qx
             zr zr or zr Rx  Ry
             zs zs zs os Sy -Sx];
    [c, n] = general.clsq(A, 2);

    B = [n [-n(2) n(1)]'];
    X = -B* [c([1 3 3 1])'; c([2 2 4 4])'];

    sq = X';

end