function PS = pointSources2Elements_faster(x,y,mesh)
%pointSources2Elements_faster generates point-wise data for mesh.
%
% PS = pointSources2Elements_faster(x,y,mesh)
% generates interpolation data array PS, where
% PS(1,:) = x
% PS(2,:) = y
% PS(3,:) = element of mesh where each (x,y) belongs

Ne = size(mesh.t,2);
Ns = max(size(x));
PS = zeros(3, Ns);

for k = 1:Ne
    I = inpolygon(x, y, mesh.p(1, mesh.t(:,k))', mesh.p(2, mesh.t(:,k))');
    if any(I)
        I = find(I);
        PS(:, I) = [x(I); y(I); k*ones(1, numel(I))];
        %PS(:,kp) = [x(kp);y(kp);k];
    end
end


end