classdef (Sealed) ArenaShape
    % ArenaShape - Class to represent arena shape types
    %
    properties (Constant)
        NoShape = NaN
        Box = 1
        Circle = 2
        Cylinder = 2
        Track = 3

        % mind the order of elements!
        ShapesAsStrings = {'No shape', 'Box', 'Circle/Cylinder', 'Linear track'};
        ListOfProperties = {'Box', 'Circle', 'Track'};
    end

    methods(Static)
        function res = allShapes()
            % Get vector with all possible shape values
            res = [helpers.ArenaShape.Box helpers.ArenaShape.Track helpers.ArenaShape.Circle];
            res = sort(res);
        end

        function res = stringForValue(value)
            res = '';
            if isnan(value)
                res = helpers.ArenaShape.ShapesAsStrings{1};
            end
            properties = helpers.ArenaShape.ListOfProperties;
            for i = 1:numel(properties)
                if helpers.ArenaShape.(properties{i}) == value
                    res = helpers.ArenaShape.ShapesAsStrings{i+1};
                end
            end
        end

        function res = valueForString(str)
            res = nan;
            properties = helpers.ArenaShape.ListOfProperties;
            idx = find(strcmpi(helpers.ArenaShape.ShapesAsStrings, str));
            if isempty(idx) || idx == 1
                return;
            end
            res = helpers.ArenaShape.(properties{idx-1});
        end

        %ShapesAsValues = [helpers.ArenaShape.NoShape helpers.ArenaShape.Box helpers.ArenaShape.Track helpers.ArenaShape.Circle];
    end
end
