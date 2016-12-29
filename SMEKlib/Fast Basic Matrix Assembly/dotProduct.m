function E = dotProduct(V1, V2)
%dotProduct vectorized dot product.
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University


E = sum(V1.*V2, 1);
%E = sum(bsxfun(@times, V1, V2), 1);
end