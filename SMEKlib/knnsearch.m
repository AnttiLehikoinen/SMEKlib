function I = knnsearch(X1, X2)

I = zeros( size(X2,1), 1);

for k = 1:size(X2,1)
    dists = sum(bsxfun(@minus, X1, X2(k,:)).^2, 2);
    [~, ind] = min(dists);
    I(k) = ind;
end