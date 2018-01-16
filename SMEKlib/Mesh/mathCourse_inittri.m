%----------------------------------------------------------------- 
%
% Antti Hannukainen and Mika Juntunen 3.11.2005
%
%-----------------------------------------------------------------
%
% Construct a mesh - structure from a given triangulation (p,t
% explained in detail alongside the mesh structure) 
% 
% TRIMESH STRUCTURE :
% 
%   p       = nodes in a 2xN-matrix
%
%   t       = triangles in a 3xN- matrix (nodes of the triangles as indeces to the p- matrix).
% 
%   edges   = a matrix of all edges in the mesh. Each column is an edge :
%             [n1 ; n2] where n1 < n2;
%
%   t2e     = a matrix connecting  triangle's and edges's. 
%             Each column corresponds to a triangle and has
%             triangle's edges in the order n1->n2, n2->n3,
%             n1->n3. 
%
%   e2t     = inverse of t2e.
%
%
% CALLING SYNTAX IS 
%
% function mesh = inittri(p,t)
%
%   p       = nodes
%   t       = triangles
%     
%   mesh    = trimesh structure corresponding to (p,t)
%


function mesh = mathCourse_inittri(p,t)

[t, ] = sort(t(1:3,:));

mesh.p  = p;
mesh.t  = t;

% Initalize size variables
Nt = size(t,2);

e = [1 2; 2 3; 1 3]';

edges = [];
for i=1:size(e,2)
  edges = [edges [ sort( [t(e(1,i),:); t(e(2,i),:)],1)]];
end

[mesh.edges,~,mesh.t2e] = sc_unique(edges);
mesh.t2e = reshape(mesh.t2e,Nt,3)';

% mesh.e2t
e = [mesh.t2e(1,:) mesh.t2e(2,:)  mesh.t2e(3,:)];
t = repmat([1:Nt],1,3);

[ef,If]= unique(e,'first');
[el,Il]= unique(e,'last');

mesh.e2t(1,ef) = t(If);
mesh.e2t(2,el) = t(Il);

mesh.e2t(2, (mesh.e2t(1,:)-mesh.e2t(2,:))==0)=0;

end



function [B,I,J] = sc_unique(A)

% sort columnwise
A = sort(A,1);

% unique
[B,I,J] = unique(A','rows');

% transpose to columns + "unsort"
B = B';

I = I(:)';
J = J(:)';

end