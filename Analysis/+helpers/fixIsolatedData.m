% Fix data with NaN values in first/last samples to improve interpolation results
%
% If there are undefined points right in the beginning of data (or in the end),
% then interpolation will be very poor. This function will assign the first
% non-NaN value (minus some offset) to the first data point. This helps
% interpolation.
function positions = fixIsolatedData(positions)
    numColumns = size(positions, 2);
    for c = 1:numColumns
        numPoints = length(positions(:, c));
        
        if isnan(positions(1, c))
            j = 2;
            while (j < numPoints) && isnan(positions(j, c))
                j = j + 1;
            end
            % do not scale if j is too big, this will create artefacts
            if j < 1000
                positions(1, c) = positions(j, c) - (0.02 * (j - 1)); % scale to the number of missed points.
                                                            % 0.02 is just a number without any particular meaning.
            else
                positions(1, c) = positions(j, c);
            end
        end
        
        if isnan(positions(numPoints, c))
            j = 2;
            while ((numPoints-j) > 1) && isnan(positions(numPoints - j, c))
                j = j + 1;
            end
            if j < 1000
                positions(numPoints, c) = positions(numPoints - j, c) + (0.02 * (j - 1));
            else
                positions(numPoints, c) = positions(numPoints - j, c);
            end
        end
    end
end
