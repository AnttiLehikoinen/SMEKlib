%----------------------------------------------------------------- 
%
% Antti Hannukainen 30.8.2007 Oulunkylä
%
%-----------------------------------------------------------------
%
% 
% Uniformly refine a triangular mesh.
%
% CALLING SYNTAX IS 
%
% function rmesh = refinetri(mesh)
% 
%

function rmesh = mathCourse_refinetri(mesh)

t = mesh.t;
p = mesh.p;
e = mesh.edges;
t2e = mesh.t2e;

Nt = size(t,2);
Ne = size(e,2);

for n=1:2
  e_nodes(n,:) = sum( [p(n,e(1,:)) ; p(n,e(2,:))])/2;
end

% Create new mesh
rmesh.p = [p  e_nodes];


% Create New elements
eb = size(p,2); 

% Edges as n1->n2, n2->n3, n1->n3.
new_t = [t(1,:) ; t2e(1,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [new_t];

new_t = [t(2,:) ; t2e(1,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t(3,:) ; t2e(3,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t2e(1,:)+eb ; t2e(2,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [rmesh.t  new_t];

rmesh = mathCourse_inittri(rmesh.p,rmesh.t);
end