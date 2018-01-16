function x = toColumn(x)
%toColumn guess what.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

if ~iscolumn(x)
    x = transpose(x);
end

end