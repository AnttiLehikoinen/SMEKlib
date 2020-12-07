function I = knnsearch(X1, X2, varargin)

I = zeros( size(X2,1), 1);

if numel(varargin)
    weights = varargin{1};
    X1 = X1 .* weights;
    X2 = X2 .* weights;
end

for k = 1:size(X2,1)
    dists = sum(bsxfun(@minus, X1, X2(k,:)).^2, 2);
    [~, ind] = min(dists);
    I(k) = ind;
end