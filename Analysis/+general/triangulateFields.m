function dt = triangulateFields(fields)
    x = [fields(:).x]';
    y = [fields(:).y]';
    dt = delaunayTriangulation([x y]);
end