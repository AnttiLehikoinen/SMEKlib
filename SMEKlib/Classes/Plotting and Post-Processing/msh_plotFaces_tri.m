function [] = msh_plotFaces_tri(msh, faces, varargin)

%generating face
if faces < 0
    faces = 1:size(msh.faces_tri,2);
end

edges = toRow(unique(abs( msh.faces_tri(:, faces) )));
msh_plotEdges3D(msh, edges, varargin{:});

%{
n = zeros(3, size(faces,2));
for kn = 1:3
    n(kn, faces(kn,:)>0) = msh.edges(1, faces(kn, faces(kn,:)>0));
    n(kn, faces(kn,:)<0) = msh.edges(2, -faces(kn, faces(kn,:)<0));
end
n = [n; n(1,:)];

plot3( reshape(msh.nodes(1, n), 4, []), ...
    reshape(msh.nodes(2, n), 4, []), ...
    reshape(msh.nodes(3, n), 4, []), ...
    varargin{:});
%}

end