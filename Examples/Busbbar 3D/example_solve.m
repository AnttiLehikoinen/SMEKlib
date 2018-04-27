%Solve busbar example.
% 
% (c) 2018 Antti Lehikoinen / Aalto University

%generating edges
edgeDefs = [1 2; 1 3; 1 4; 2 3; 3 4; 4 2]; %numbering of edges in reference element

edges = [t(edgeDefs(1,:),:) t(edgeDefs(2,:),:) t(edgeDefs(3,:),:) t(edgeDefs(4,:),:) t(edgeDefs(5,:),:) t(edgeDefs(6,:),:)];
edges = sort(edges,1);
edges = unique(edges', 'rows')';

%elements to edges array
Ne = size(t,2);
t2e = zeros(6, Ne);
for k = 1:6
    edges_el = t(edgeDefs(k,:), :);
    [~, loc] = ismember(sort(edges_el,1)', edges', 'rows');
    t2e(k, :) = loc;
    
    inds_rev = find( edges_el(1,:) ~= edges(1,loc) );
    t2e(k, inds_rev) = -1*t2e(k, inds_rev);
end

%getting free variables
t_dir = Surfaces.get('Outer');
e_dir = [t_dir([1 2],:) t_dir([2 3],:) t_dir([3 1],:)];
[~, loc] = ismember(sort(e_dir,1)', edges', 'rows');
e_dir = loc';
e_free = setdiff(1:size(edges,2), e_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mesh
msh = SimpleTetMesh(p, t, edges, t2e);

Wcurl = Nedelec3D(Operators.curl); %shape function
W = Nedelec3D(Operators.I);

%assembling brutely-invertible matrix
Mc = MatrixConstructor();
Mc.assemble_matrix(Wcurl, Wcurl, 1, [], msh);
Mc.assemble_matrix(W, W, 1e-9, [], msh); %uglyyyy

S = Mc.finalize();

%current source
Fc = MatrixConstructor();
els = Surfaces.get('Bar');
J = 1000*repmat( [0;0;1], 1, numel(els));

F = Fc.assemble_vector(W, 1,  J, els, msh).finalize();

figure(2); clf;
spy(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct solution, will overfill RAM on larger problems.
A = zeros(size(msh.e,2), 1);

tic
%A(e_free) = S(e_free,e_free) \ F(e_free);
toc

%plotting solution at element centers
Avec = zeros(3, Ne);
xref = [0.25; 0.25; 0.25];
for k = 1:6
    Avec = Avec + bsxfun(@times, Wcurl.eval(k, xref, msh, []), A(abs(msh.t2e(k,:)))' );
    %Avec = Avec + Wcurl.eval(k, xref, msh, []);
end
Avec = real(Avec);

Xplot = zeros(3, Ne);
for k = 1:4
    Xplot = Xplot + msh.p(:, msh.t(k,:));
end
Xplot = Xplot/4;

figure(3); clf; hold on;
quiver3(Xplot(1,:), Xplot(2,:), Xplot(3,:), Avec(1,:), Avec(2,:), Avec(3,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterative solution

%assembling the (singular) matrix. Seems to work fine here, but in general
%the load vector F should be projected to the rage of S, i.e. making sure
%it's divergence free.
Mc = MatrixConstructor();
Mc.assemble_matrix(Wcurl, Wcurl, 1, [], msh);

S = Mc.finalize();

%preconditioning
tic;
%[L,U] = ilu( S(e_free,e_free), struct('type', 'ilutp', 'droptol', 0.01) );
[L,U] = ilu( S(e_free,e_free) );
toc

%solving
%x = bicgstabl(S(e_free,e_free), F(e_free), 1e-6, 100, L, U, []);
x = gmres(S(e_free,e_free), F(e_free) , 20, 1e-6, 10, L, U);
toc
A2 = zeros(size(A));
A2(e_free) = real(x);

%plotting
Avec2 = zeros(size(Avec));
for k = 1:6
    Avec2 = Avec2 + bsxfun(@times, Wcurl.eval(k, xref, msh, []), A2(abs(msh.t2e(k,:)))' );
    %Avec = Avec + Wcurl.eval(k, xref, msh, []);
end

figure(4); clf; hold on;
quiver3(Xplot(1,:), Xplot(2,:), Xplot(3,:), Avec2(1,:), Avec2(2,:), Avec2(3,:));