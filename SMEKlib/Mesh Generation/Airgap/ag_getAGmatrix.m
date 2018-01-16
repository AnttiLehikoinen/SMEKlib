function S = ag_getAGmatrix(this, rotorAngle, varargin)

if numel(varargin)
    N = varargin{1};
else
    N = size(this.S_const, 1);
end
Np = size(this.S_const, 1);

%...stuff goes here

%rotating virtual nodes
%{
rm = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)];
this.msh_ag.p(:, this.n_moving) = rm*this.p_virt(:, this.n_moving);

%shifting indices in triangulation
nodeShift = -floor( (rotorAngle - this.shiftTol/2) / this.shiftTol ) - 1;
newPositions = mod(this.original_positions - 1 + nodeShift, numel(this.n_bnd) ) + 1;

t_moving = this.t_moving;
t_moving(this.inds_r) = this.n_bnd(newPositions);

this.msh_ag.t = t_moving;
%}

[t_moving, p_virt] = this.update(rotorAngle, varargin{:});

this.msh_ag.t = t_moving;
this.msh_ag.p = p_virt;

Sag_c = MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), 1/(pi*4e-7), [], this.msh_ag);

%moving to global indexing and taking care of symmetry sectors
inds = 1:(Sag_c.Nvals);
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.I(inds));
Sag_c.I(inds) = this.el_table(2, Sag_c.I(inds));
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.J(inds));
Sag_c.J(inds) = this.el_table(2, Sag_c.J(inds));

%finalizing
S = blkdiag(this.S_const, sparse(N-Np, N-Np)) + Sag_c.finalize(N, N);
end