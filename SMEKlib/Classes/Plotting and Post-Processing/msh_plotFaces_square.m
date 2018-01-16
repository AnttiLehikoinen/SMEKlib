function [] = msh_plotFaces_square(msh, faces, varargin)

%generating face
if faces < 0
    faces = 1:size(msh.faces_square,2);
end

edges = toRow(unique(abs( msh.faces_square(:, faces) )));
msh_plotEdges3D(msh, edges, varargin{:});

end