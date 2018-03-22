function this = ag_setConstantMatrix(this, Np)
%ag_setConsttantMatrix set constant part of airgap matrix.
%
% (c) 2018 Antti Lehikoinen / Aalto University

%generating constant-part of ag matrix
msh = SimpleMesh(this.p_virt, this.t_const);

Sag_c = ...
    MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), 1/(pi*4e-7), [], msh);

%moving to global indexing and taking care of symmetry sectors
inds = 1:Sag_c.Nvals;
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.I(inds));
Sag_c.I(inds) = this.el_table(2, Sag_c.I(inds));
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.J(inds));
Sag_c.J(inds) = this.el_table(2, Sag_c.J(inds));
this.S_const = Sag_c.finalize(Np,Np);

end
