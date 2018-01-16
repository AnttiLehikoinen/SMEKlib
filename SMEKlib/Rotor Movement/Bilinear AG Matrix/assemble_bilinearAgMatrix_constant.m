function Sag = assemble_bilinearAgMatrix_constant(msh, dims)
%assemble_bilinearAgMatrix_constant assembles the constant components of
%the bilinear air-gap stiffness matrix.
% 
% Entries
% Sag_ij = Int{ B_i x B_j r*dtheta dr,
% with
% B_i = phi_i(theta) * p(r)
% where phi_i = trace of the shape function N_i on the air-gap boundary,
% and p(r) is a 1st-order polynomial satisfying p(rr)=1, pr(rs)=0 if the
% node i belongs to the rotor, and p(rr)=0, p(rs)=1 otherwise.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

N_ag_s = msh.bndData.N_ag_s;
N_ag_r = msh.bndData.N_ag_r;
N_ag = N_ag_s + N_ag_r;

agEdges = [1:N_ag_s (N_ag_s+1):N_ag;
    2:N_ag_s 1 (N_ag_s+2):N_ag (N_ag_s+1)];

%signed edge lengths in radians
le = angleDifference( msh.bndData.agAngles_global(agEdges(2,:)), msh.bndData.agAngles_global(agEdges(1,:)) );
Ne = size(agEdges, 2);

%reference 1D shape functions
pref = { [0 1]; [1 -1] };
dpref = { [1 0]; [-1 0] };
NO = 2;

ri = dims.D_ro / 2;
ro = dims.D_si / 2;

%assembling indices and entries
[x_quad, w_quad] = get_1DQuadPoints(4); %Gaussian quadrature points and weights
x_quad = (1+x_quad)/2; w_quad = w_quad / 2; %scaling to the reference interval [0,1]

N_quad = numel(w_quad);
I = zeros(1, Ne*NO^2); J = zeros(1, Ne*NO^2); E = zeros(1, Ne*NO^2);

% calculating entries

%integrating radial function terms (constant over circumference)
%lazy implementation for now in global coordinates
radii = linspace(ri, ro);
pr = 1 - (radii-ri)./(ro-ri);
ps = (radii-ri)./(ro-ri);

%p_dr term (same for both stator and rotor
p_dr = repmat(1/(ro-ri)^2 * trapz(radii, radii), 1, Ne);
%p_dtheta term
p_dtheta = [repmat(trapz(radii, ps.^2./radii), 1, N_ag_s) repmat(trapz(radii, pr.^2./radii), 1, N_ag_r)];

for k_quad = 1:N_quad
    ri = 1;
    p_quad = [1; x_quad(k_quad)];
    for k_test = 1:NO
        for k_shape = 1:NO
            inds = ((ri-1)*Ne + 1):(ri*Ne);
            
            %dr term
            E(inds) = E(inds) + ...
                w_quad(k_quad)*(pref{k_test}*p_quad)*(pref{k_shape}*p_quad) .* p_dr .* abs(le);
            
            %dtheta term
            E(inds) = E(inds) + ...
                w_quad(k_quad)*(dpref{k_test}*p_quad./le).*(dpref{k_shape}*p_quad./le) .* p_dtheta .* abs(le);
                
            
            %setting indices
            I(inds) = agEdges(k_test,:);
            J(inds) = agEdges(k_shape,:);
            
            ri = ri + 1;
        end
    end
end
%cleaning up (not sure if necessary)
I = I(E~=0); J = J(E~=0); E = E(E~=0);

%handling periodicity conditions if necessary
if isfield(msh, 'symmetrySectors')      
    E = E .* (msh.periodicityCoeff.^msh.bndData.el_table(3, I));
    I = msh.bndData.el_table(2, I);
    E = E .* (msh.periodicityCoeff.^msh.bndData.el_table(3, J));
    J = msh.bndData.el_table(2, J);
    
    E = E / msh.symmetrySectors;
    
    N_ag_true = msh.bndData.N_ag_s_true + msh.bndData.N_ag_r_true;
else
    N_ag_true = N_ag;
end

Sag = sparse(I, J, E, N_ag_true, N_ag_true);

end