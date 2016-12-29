function detF = mappingDeterminant(F)
%mappingDeterminant determinant of the affine mapping.
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

detF = F(1,:).*F(4,:) - F(2,:).*F(3,:);
end