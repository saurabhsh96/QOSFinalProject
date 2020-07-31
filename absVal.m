%Finding magnitude
function [mag, maxim] = absVal(components)
    mag = sqrt((abs(components{1}).^2) + (abs(components{2}).^2) ...
        + (abs(components{3}).^2));
    maxim = max(max(mag));
end