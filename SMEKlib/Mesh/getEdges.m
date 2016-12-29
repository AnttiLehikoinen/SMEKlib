function [edges, e2t, t2e] = getEdges(t)
%getEdges returns the edge definition.
%
% [edges, e2t, t2e] = getEdges(t)
% when given a 3xNe triangulation t, the function returns the following
% arrays
%   edges = 2xNedges array, containing the start- and end-nodes of each
%      edge in the mesh (on first and second row respectively)
%   e2t = a 2xNedges array, containing the triangles that the mesh belongs
%       to
%   t2e = a 3xNe array, containing the edges that each triangle has
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University,
%   based on the work of Antti Hannukainen and Mika Juntunen, 3.11.2005

%

% initialization
t = sort(t,1);
Nt = size(t, 2);

%edges
edges = sort([t(1,:) t(2,:) t(3,:);
    t(2,:) t(3,:) t(1,:)], 1);

[edges, ~, t2e] = sc_unique(edges);
t2e = reshape(t2e,Nt,3)';

%getting the edges-to-elements array
e = [t2e(1,:) t2e(2,:)  t2e(3,:)];
t_list = repmat(1:Nt, 1,3);

[ef,If]= unique(e, 'first');
[el,Il]= unique(e, 'last');

e2t(1,ef) = t_list(If);
e2t(2,el) = t_list(Il);

e2t(2, (e2t(1,:)-e2t(2,:))==0)=0;

end



function [B,I,J] = sc_unique(A)

% sort columnwise
A = sort(A, 1);

% unique
[B,I,J] = unique(A','rows');

%untranspose
B = B';

%unrolling
I = I(:)';
J = J(:)';

end