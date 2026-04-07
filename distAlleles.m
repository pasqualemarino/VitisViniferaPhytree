function A = distAlleles(S, distanceName) %S is matrix without header
    if iscell(S)
        S = cell2mat(S);
    end
    S = double(S);
    A = pdist2(S', S', distanceName);

end

