% Determine sample time from position data
%
function sampleTime = sampleTimeFromData(positions)
    i = 3;
    while isnan(positions(i, 1)) || isnan(positions(i-1, 1))
        i = i + 1;
    end
    sampleTime = positions(i, 1) - positions(i-1, 1);
end
