function val = aux_polyvals(ps, xs)
%val = aux_polyvals(ps, xs) returns the values val of polynomials ps evaluated at the points xs
%   val = N x n array of polynomial values
%   ps = N x (p+1) array of polynomial coefficients, of order p
%   xs = an array of size either N x n or 1 x n, in the latter case all
%   polynomials are evaluated at the same point

N = size(ps, 1);
n = size(xs, 2);
pp = size(ps, 2);

val = zeros(N, n); xtemp = ones(size(xs));
for k_term = 1:pp
    %val = val + bsxfun(@times, ps(:,k_term), xs.^(pp - k_term));
    
    %val = val + bsxfun(@times, ps(:,k_term), realpow(xs, (pp - k_term)));
    
    val = val + bsxfun(@times, ps(:,pp - k_term + 1), xtemp); 
    %val = val + aux_bsxtimes(ps(:,pp - k_term + 1), xtemp); 
    xtemp = xtemp .* xs;
    
    %temp = xs.^(pp - k_term);
    %val = val + aux_bsxtimes(ps(:,k_term), temp);
end

end

function y = aux_bsxtimes(x1, x2)

y = zeros(size(x1,1), size(x2,2));
for col = 1:size(x2,2)
    y(:,col) = x1.*x2(:,col);
end


end