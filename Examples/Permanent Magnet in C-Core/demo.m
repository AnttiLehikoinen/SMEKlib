%simulates a simple PM in C-core problem
%
% (c) 2018 Antti Lehikoinen / Aalto University


addpath(genpath('..\..\SMEKlib'));
gmsh_path = 'E:\Software\Work\gmsh'; %CHANGE this


%loading problem
geo_path = [pwd '\PM.geo'];
gw = gwrap(gmsh_path); %creating gwrap object
[p, t, Surfaces] = gw.loadMesh(geo_path); %loading mesh

msh = SimpleMesh(p,t);

%plotting
figure(1); clf; hold on; box on; axis equal;
msh_fill(msh, Surfaces.get('Core'), [1 1 1]*0.7);
msh_fill(msh, Surfaces.get('Magnet'), 'r');
msh_fill(msh, Surfaces.get('Air'), 'g');

edges_out = Surfaces.get("Out");
n_Dirichlet = toRow(unique(edges_out));
msh_plot(msh, n_Dirichlet, 'ko');

%solving
Np = size(msh.p,2);
nfree = setdiff(1:Np, n_Dirichlet);
nu0 = 1/(pi*4e-7);
nu = nu0*ones(1, size(msh.t,2)); %reluctivity vector
nu(Surfaces.get('Core')) = nu0/5000;

Ncurl = Nodal2D(Operators.curl); %curl of nodal shape function
S = MatrixConstructor(Ncurl, Ncurl, nu, [], msh).finalize();

%load vector from PMs
el_pm = Surfaces.get('Magnet');
Br = 1.4; %permanent magnet
Hrem = 1/(pi*4e-7)*repmat([0; Br], 1, numel(el_pm));
F_PM = MatrixConstructor().assemble_vector(Ncurl, 1, Hrem, el_pm, msh).finalize();

A = zeros(Np, 1);
A(nfree) = S(nfree,nfree) \ F_PM(nfree);

figure(2); clf; hold on; box on;
drawFluxDensity(msh, A, 'LineStyle', 'none'); 
colormap('jet'); colorbar; caxis([0 2])
drawFluxLines(msh, A, 16, 'k');
axis equal;