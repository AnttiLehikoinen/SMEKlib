function [] = msh_plot3D(msh, nodes, varargin)

plot3(msh.nodes(1,nodes), msh.nodes(2,nodes), msh.nodes(3,nodes), varargin{:});

end