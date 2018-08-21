%Example of eddy currents in a steel sheet.
%
% Loads the mesh and solves a linear time-harmonic problem.
%
% (c) 2018 Antti Lehikoinen / Aalto University

addpath(genpath('..\..\SMEKlib'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load mesh and simulate a time-harmonic example
gw = gwrap('dummy_not_needed');
[p,t,Surfaces] = gw.loadMesh([cd() '\sheet.msh']);


msh = MachineMesh(p, t);

%ordering edges as Nedelec elements need them like that
edgeDefs = [2 3; 3 1; 1 2];

%re-definition of edges to comply with the definition of Nedelec elements,
%should be fixed
edges = [msh.t(edgeDefs(1,:),:) msh.t(edgeDefs(2,:),:) msh.t(edgeDefs(3,:),:)];
edges = sort(edges,1);
edges = unique(edges', 'rows')';
msh.edges = edges;

%getting elements-to-edges array with directions indicated by the sign
Ne = size(t,2);
t2e = zeros(3, Ne);
for k = 1:3
    edges_el = msh.t(edgeDefs(k,:), :);
    [~, loc] = ismember(sort(edges_el,1)', msh.edges', 'rows');
    t2e(k, :) = loc;
    
    inds_rev = find( edges_el(1,:) ~= msh.edges(1,loc) );
    t2e(k, inds_rev) = -1*t2e(k, inds_rev);
end
msh.t2e = t2e;

%finding outer edges
e_dir = Surfaces.get('Out');
[~, loc] = ismember(sort(e_dir,1)', msh.edges', 'rows');
e_dir = loc';
%e_dir = e_dir(1);
e_free = setdiff(1:size(msh.edges,2), e_dir);

e_free = setdiff(1:size(msh.edges,2), []); %keep all free for laziness

figure(1); clf; hold on; axis equal;
msh_triplot(msh, [], 'b');
msh_triplot(msh, [154 448], 'r');
%plot( [msh.p(1,e_dir(1,:)); msh.p(1,e_dir(2,:))], [msh.p(2,e_dir(1,:)); msh.p(2,e_dir(2,:))], 'r');
plot( [msh.p(1,msh.edges(1,e_dir)); msh.p(1,msh.edges(2,e_dir))], [msh.p(2,msh.edges(1,e_dir)); msh.p(2,msh.edges(2,e_dir))], 'r');

ke = 17;
plot( [msh.p(1,msh.edges(1,ke)); msh.p(1,msh.edges(2,ke))], [msh.p(2,msh.edges(1,ke)); msh.p(2,msh.edges(2,ke))], 'ko-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing
Phi = 1.5*1e-3*15e-3; %total flux, corresponding to average flux density
nu = 1/(5000 * pi*4e-7); %reluctivity
sigma = 1e6; %conductivity
w = 2*pi*500; %angular frequency

%assembling system
Wcurl = Nedelec2D(Operators.curl); %Nedelec-curl function
W = Nedelec2D(Operators.I); %Nedelec function
S_AA = MatrixConstructor(Wcurl, Wcurl, nu, [], msh).finalize(); %stiffness matrix
M_AA = MatrixConstructor(W, W, sigma, [], msh).finalize(); %mass matrix

%constraint vector for total flux
V_Constr = MatrixConstructor().assemble_vector(Wcurl, 1, 1, [], msh).finalize();

%solving
Q = [S_AA(e_free, e_free)+1i*w*M_AA(e_free,e_free) V_Constr(e_free);
    V_Constr(e_free)' 0];
X = Q \ [zeros(numel(e_free),1); Phi];

A = zeros(size(msh.edges, 2),1);
A(e_free) = X(1:numel(e_free));

%plotting at element centers
Jvec = zeros(2, Ne);
Avec = zeros(1, Ne);
xref = [0.25; 0.25];
for k = 1:3
    Jvec = Jvec + bsxfun(@times, W.eval(k, xref, msh, []), 1i*sigma*w*A(abs(msh.t2e(k,:)))' );
    Avec = Avec + bsxfun(@times, Wcurl.eval(k, xref, msh, []), A(abs(msh.t2e(k,:)))' );
end
Jvec = real(Jvec);

Xplot = zeros(2, Ne);
Xtri = zeros(3, Ne); Ytri = zeros(3, Ne);
for k = 1:3
    Xplot = Xplot + msh.p(:, msh.t(k,:));
    Xtri(k,:) = msh.p(1, msh.t(k,:));
    Ytri(k,:) = msh.p(2, msh.t(k,:));
end
Xplot = Xplot/3;

%plotting flux density
figure(3); clf; hold on; axis equal tight;
fill(Xtri, Ytri, real(Avec), 'linestyle', 'none'); 
colormap('jet'); colorbar; caxis([0 1.5])

%plotting eddy currents
figure(4); clf; hold on; axis equal;
quiver(Xplot(1,:), Xplot(2,:), Jvec(1,:), Jvec(2,:));
