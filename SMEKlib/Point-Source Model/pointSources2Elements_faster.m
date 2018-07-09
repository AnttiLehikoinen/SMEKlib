function [points, elements] = pointSources2Elements_faster(x,y,mesh)
%pointSources2Elements_faster generates point-wise data for mesh.
%
% PS = pointSources2Elements_faster(x,y,mesh)
% generates interpolation data array PS, where
% points(1,:) = x
% points(2,:) = y
% elements = element of mesh where each (x,y) belongs

Ne = size(mesh.t,2);
Ns = max(size(x));
points = zeros(2, Ns);
elements = zeros(1, Ns);

for k = 1:Ne
    I = inpolygon(x, y, mesh.p(1, mesh.t(:,k))', mesh.p(2, mesh.t(:,k))');
    if any(I)
        I = find(I);
        points(:, I) = [x(I); y(I)];
        elements(:,I) =  k*ones(1, numel(I));
        %PS(:,kp) = [x(kp);y(kp);k];
    end
end


end