function [] = msh_plotFace_tri(msh, faces, varargin)

%generating face
n = zeros(3, size(faces,2));
for kn = 1:3
    n(kn, faces(kn,:)>0) = msh.edges(1, faces(kn, faces(kn,:)>0));
    n(kn, faces(kn,:)<0) = msh.edges(2, -faces(kn, faces(kn,:)<0));
end

n

plot3( reshape(msh.nodes(1, n), 1, []), ...
    reshape(msh.nodes(2, n), 1, []), ...
    reshape(msh.nodes(3, n), 1, []), ...
    varargin{:});

end