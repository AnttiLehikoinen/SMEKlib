function y = aux_logspace(x1, x2, N)
% y = aux_logspace(x1, x2, N)
% returns N points spaced logarithmically on the interval
% [x1 x2]

lob = logspace( 0, 1, N); %basic log space

y = bsxfun(@plus, bsxfun(@times, (x2-x1)/9, lob), x1 - (x2-x1)/9);

end