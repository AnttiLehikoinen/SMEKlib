function v = getShortAxis(msh, els)
%getShortAxis Get the short principal axis of the domain consisting of the
%elements "els".
%
% Based on the work
%   https://www.anttilehikoinen.fi/research-work/determining-direction-rectangle/
%
% (c) 2018 Antti Lehikoinen / Aalto University

%load vector and mass matrix
f = MatrixConstructor().assemble_vector(Nodal2D(Operators.I), 1, 1, els, msh).finalize();
M = MatrixConstructor(Nodal2D(Operators.I), Nodal2D(Operators.I), 1, els, msh).finalize();

%removing mean
x0 = msh.p*f / (ones(1, size(msh.p,2))*f);
xc = bsxfun(@plus, msh.p, -x0);

%inertia matrix for the specified domain
Mprime = xc * M * xc';

%getting principal axes
[V,D] = eig(Mprime);

%returning shorter axis
[~, ind] = min( diag(D) );
v = V(:,ind);

%v = 1;

end