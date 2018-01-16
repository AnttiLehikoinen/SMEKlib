function MF = assemble_AGbndFourierCoefficientMatrix(Nterms, msh)
% MF = assemble_AGbndFourierCoefficientMatrix(Nterms, bndData) returns a
% matrix MF for calculating the Fourier coefficients of the vector
% potential on the stator and rotor surface, w.r.t. the angular coordinate.
%
% MF has the entries
% MF_{r,c} = Int_surface{ N_c f_r(theta) dtheta }

bndData = msh.bndData;

N_ag_s = bndData.N_ag_s;
N_ag_r = bndData.N_ag_r;
N_ag = N_ag_s + N_ag_r;

edges = [1:N_ag_s (N_ag_s+1):N_ag;
    2:N_ag_s 1 (N_ag_s+2):N_ag N_ag_s+1];

Nedges = size(edges,2);

%determining shape function coefficients in the global angular coordinates
ae_starts = bndData.agAngles_global(edges(1,:))';
ae_ends = bndData.agAngles_global(edges(2,:))';

theta_0 = ae_starts;
ae_ends = angleDifference(ae_ends, ae_starts);
ae_starts = zeros(size(ae_starts));


fun_coeffs = {[1+ae_starts./(ae_ends-ae_starts) -1./(ae_ends-ae_starts)] ...
    [-ae_starts./(ae_ends-ae_starts) 1./(ae_ends-ae_starts)]};

%assembling
E = zeros(2*Nedges*Nterms, 1); I = E; J = E;

[mints_1, mints_x] = internal_integrateMonomials(ae_starts, ae_ends, theta_0, Nterms);

lin_inds = 1:Nedges*Nterms;
row_inds = mod( lin_inds-1, Nedges) + 1;
col_inds = floor( (lin_inds-1)/Nedges) + 1 + (row_inds > N_ag_s)*Nterms;
for kfun = 1:2
    inds = ((kfun-1)*Nedges*Nterms + 1):(kfun*Nedges*Nterms);
    
    Etemp = bsxfun(@times, mints_1, fun_coeffs{kfun}(:,1)) + ...
        bsxfun(@times, mints_x, fun_coeffs{kfun}(:,2));
    
    %Etemp = mints_1;
    
    E(inds) = Etemp(:);
    I(inds) = edges(kfun, row_inds);
    J(inds) = col_inds;
end

I = I(E~=0); J = J(E~=0); E = E(E~=0);

if isfield(msh, 'symmetrySectors')
    E = E .* transpose(msh.periodicityCoeff.^msh.bndData.el_table(3, I));
    I = msh.bndData.el_table(2, I);
    MF = sparse(J, I, E/pi, 2*Nterms, numel(msh.bndData.agNodes_global));
else
    MF = sparse(J, I, E/pi, 2*Nterms, N_ag);
end


end

function [mints_1, mints_x] = internal_integrateMonomials(starts, ends, theta_0, Nterms)

ns = floor( (1:Nterms) / 2 );

starts = internal_2column(starts);
ends = internal_2column(ends);

thetas = repmat(mod(1:Nterms, 2)*pi/2, numel(starts), 1);
xstarts = bsxfun(@times, starts + internal_2column(theta_0), ns);
xends = bsxfun(@times, ends + internal_2column(theta_0), ns);


mints_1 = -(cos(xends + thetas) - cos(xstarts + thetas));
mints_1 = bsxfun(@times, mints_1, 1./ns);

mints_x = bsxfun(@times, repmat(ends,1,Nterms).*cos(xends+thetas) - repmat(starts,1,Nterms).*cos(xstarts+thetas), -1./ns) + ...
    bsxfun(@times, sin(xends+thetas) - sin(xstarts+thetas), 1./ns.^2);

mints_1(:,1) = internal_2column(ends - starts);
mints_x(:,1) = 0.5*internal_2column(ends.^2 - starts.^2);

%mints_1 = sin(repmat(ends,1,Nterms) + thetas);

end

function v = internal_2column(v)

if ~iscolumn(v)
    v = transpose(v);
end

end