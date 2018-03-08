function [Force, Torque] = compute_WST2(A, msh, n_inner, n_outer, els)
%compute_WST2 force and torque computation.
% 
% (c) 2018 Antti Lehikoinen / Smeklab

mu0 = pi*4e-7;

Ne = numel(el_ag);
Np = size(msh.p, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving weighting function

% assembling stiffness matrix
%static mesh part
els_msh = els( els <= Ne ); %elements in the static mesh
msh_comp = SimpleMesh(msh.p, msh.t(:,els_msh));
S_msh = MatrixConstructor(Noda2D(Operators.grad), Noda2D(Operators.grad), 1/mu0, els_msh, msh).finalize(Np,Np);

%airgap mesh part, if any
els_ag = els( els > Ne ) - Ne; %elements in the airgap mesh
if any(els_ag)
    [tag, pag] = msh.bandData.t_ag(0);
    msh_ag = SimpleMesh(pag, tag(:,els_ag));
    Sc_ag = MatrixConstructor(Noda2D(Operators.grad), Noda2D(Operators.grad), 1/mu0, [], msh_ag);

    %fixing indexing
    inds = 1:(Sc_ag.Nvals);
    Sc_ag.E(inds) = Sc_ag.E(inds) .* msh.bandData.el_table(3, Sc_ag.I(inds));
    Sc_ag.I(inds) = this.el_table(2, Sc_ag.I(inds));
    Sc_ag.E(inds) = Sc_ag.E(inds) .* msh.bandData.el_table(3, Sc_ag.J(inds));
    Sc_ag.J(inds) = this.el_table(2, Sc_ag.J(inds));
    
    n_ag = toRow(msh.bandData.el_table(2, tag(els_ag)));
    
    S = S_msh + Sc_ag.finalize(Np,Np);
else
    msh_ag = SimpleMesh(zeros(2,0), zeros(size(msh.t,1), 0));
    S = S_msh;
    n_ag = [];
end

%actual solution
g = zeros(Np, 1);
g(n_outer) = 1;
n_free = setdiff([reshape(msh.t(:, els_msh),1,[]) n_ag], [n_outer n_inner]);
g(n_free) = S(n_free, n_free) \ (-S(n_free, n_outer)*g(n_outer));

[x_quad, W_quad] = get_2DtriangleIntegrationPoints(1); %FIXME generalize this!
N_quad = numel(W_quad);

Torque = zeros(1, size(A,2));
Force = zeros(2, size(A,2));

N = Nodal2D(Operators.grad);
for k_quad = 1:N_quad
    %static part
    [F, F0] = msh_comp.getMappingMatrix([], x_quad(:,k_quad)); detF = mappingDeterminant(F);
    g_grad = zeros(2, size(msh_comp.t, 2));
    x_global = mappingTimesVector(x_quad(:,k_quad), false, false, F) + F0;
    
    %computing flux density
    for kn = 1:size(msh_comp.t,1)
        B = B + bsxfun(@times, N.eval(k_test, x_quad(:,k_quad), msh, F, detF), ...
            transpose(A(msh_comp.t(kn, :))) );
        g_grad = g_grad + bsxfun(@times, N.eval(k_test, x_quad(:,k_quad), msh, F, detF), ...
            transpose(g(msh_comp.t(kn, :), :)));
    end
    
    %airgap part, if any
    if els_ag
        [Fa, F0a] = msh_ag.getMappingMatrix([], x_quad(:,k_quad)); 
        detFa = mappingDeterminant(Fa);
        g_grad = zeros(2, size(msh_ag.t, 2));
        x_global = mappingTimesVector(x_quad(:,k_quad), false, false, F) + F0;
    end
    
    
    %[Fag, F0ag] = msh_ag.getMappingMatrix(els_ag, x_quad(:,k_quad));

end