function [points, inds] = gw_addTolerances(points, tol)


%array for holding the surface points
x = zeros(2, 1000); 
inds = zeros(1, 1000);
np = 0;
function addx(xp, id)
    % nested function for incrementing the array
    if (np+size(xp,2)) > size(x,2)
        nadd = min(2*size(x,2), size(xp,2));
        x = [x zeros(2,nadd)];
    end
    x(:,(np+1):(np+size(xp,2))) = xp;
    inds(:,(np+1):(np+size(xp,2))) = id;
    np = np + size(xp,2);
end

Np = size(points, 2);
for k = 1:Np
    xstart = points(:,k);
    xend = points(:, mod(k, Np)+1);
    
    Np_seg = ceil( norm(xend-xstart)/tol ) + 1;
    xp = bsxfun(@plus, xstart, bsxfun(@times, xend-xstart, linspace(0,1,Np_seg)));
    addx(xp(:,1:(end-1)), k);
end

points = x(:, 1:np);
inds = inds(1:np);

end