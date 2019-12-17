addpath(genpath('..\..\SMEKlib'));

gw = gwrap('E:/Software/Work/gmsh43');

fname = 'pmsm.geo';

%meshing
[p, t, Surfaces] = gw.loadMesh(fname);