function y = zeropadcat(varargin)

Nc = sum( cellfun('size', varargin, 2) );
Nr = max( cellfun('size', varargin, 1) );

y = zeros(Nr,Nc);
ri = 1;
for k = 1:size(varargin,2)
    [r,c] = size(varargin{k});
    y(1:r,ri:(ri+c-1)) = varargin{k};
    ri = ri+c;
end

end