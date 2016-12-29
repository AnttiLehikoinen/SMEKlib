function y = aux_linspace(x1, x2, Npoints)
%aux_linspace auxiliary linspace function
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

blin = linspace(1, 10, Npoints); 

y = bsxfun(@plus, bsxfun(@times, (x2-x1)/9, blin), x1 - (x2-x1)/9);

end