function x = toRow(x)
%toRow converts vector to row-vector.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if iscolumn(x)
    x = transpose(x);
end

end