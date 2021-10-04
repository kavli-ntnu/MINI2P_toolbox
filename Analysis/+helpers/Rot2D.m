function R = Rot2D(angle)
    R = [cos(angle), -sin(angle); 
        sin(angle), cos(angle)];
end