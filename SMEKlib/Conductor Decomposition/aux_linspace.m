function y = aux_linspace(x1, x2, Npoints)

blin = linspace(1, 10, Npoints); 

y = bsxfun(@plus, bsxfun(@times, (x2-x1)/9, blin), x1 - (x2-x1)/9);

end