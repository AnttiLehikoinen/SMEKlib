function FV = mappingTimesVector(V, inverse, transpose, F, varargin)
%mappingTimesVector applies the affine mapping on a vector.
% 
% Call syntax options
% FV = mappingTimesVector(V, inverse, transpose, F)
% 
% FV = mappingTimesVector(V, inverse, transpose, F, F0) include vector term
% F0.
% 
% Note: even with inverse=true the division by the mapping determinant is
% be default ignored; use
% FV = mappingTimesVector(V, inverse, transpose, F, F0, detF) to divide
% include it.
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

if transpose
    ind_11 = 1;
    ind_21 = 3;
    ind_12 = 2;
    ind_22 = 4;
else
    ind_11 = 1;
    ind_21 = 2;
    ind_12 = 3;
    ind_22 = 4;
end

if ~isempty(varargin) && ~isempty(varargin{1})
    V = bsxfun(@plus, V, varargin{1});
end

if inverse
    %returns vectorized det(F)*F\V or det(F)*F'\V
    FV = [F(ind_22,:).*V(1,:)-F(ind_12,:).*V(2,:);
        -F(ind_21,:).*V(1,:)+F(ind_11,:).*V(2,:)];
else
    %returns vectorized F*V or F'*B
    FV = [F(ind_11,:).*V(1,:)+F(ind_12,:).*V(2,:);
        F(ind_21,:).*V(1,:)+F(ind_22,:).*V(2,:)];
end

if numel(varargin) == 2
    FV = bsxfun(@times, FV, 1./varargin{2});
end

end