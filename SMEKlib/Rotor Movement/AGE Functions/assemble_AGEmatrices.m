function [MF, A2ag, ag2H, Q, N] = assemble_AGEmatrices(bndData, dims, Np, varargin)

% assembles the typical AGE matrices

nu0 = 1/(pi*4e-7); %air reluctivity
r_si = dims.D_si / 2;
r_ro = dims.D_ro / 2;

N_ag_s = bndData.N_ag_s;
N_ag_r = bndData.N_ag_r;

if numel(varargin)
    Nfreqs = varargin{1};
else
    Nfreqs = floor( max(N_ag_s, N_ag_r) / 2 -0.1);
end

% assembling F matrix: a matrix for computing the spatial Fourier series on
% the stator and rotor ag-surface, from the boundary potential values
MF = assemble_AGbndFourierCoefficientMatrix(2*Nfreqs + 1, bndData);
MF = MF([2:(2*Nfreqs) (2*Nfreqs+2):(4*Nfreqs+2)], :); %discarding DC components

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling A2ag matrix: a matrix for computing the air-gap solution
% coefficients from the boundary Fourier series coefficients
I = zeros(8*Nfreqs, 1); J = I; E = I;

ksi = r_si / r_ro;
for kfreq = 1:Nfreqs
    localInverse = inv([ksi^kfreq ksi^-kfreq;1 1]);
    
    inds = ((kfreq-1)*8 + 1):(kfreq*8);
    
    I(inds) = (kfreq-1)*4 + [1 2 1 2 3 4 3 4];
    J(inds) = (kfreq-1)*2 + [1*[1 1] (2*Nfreqs+1)*[1 1] 2*[1 1] (2*Nfreqs+2)*[1 1]];
    E(inds) = [localInverse(:); localInverse(:)];
end
A2ag = sparse(I, J, E, 4*Nfreqs, 4*Nfreqs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling the ag2H matrix: a matrix from computing the Fourier series
% representation of the circumferential component of H on the stator/rotor
% surface
I = zeros(8*Nfreqs, 1); J = I; E = I;

for kfreq = 1:Nfreqs
    %stator surface
    localMatrix_stator = nu0*kfreq/r_si * [ksi^kfreq*[1;-1] (ksi^-kfreq)*[1;-1]];
    inds = (kfreq-1)*8 + (1:4);
    
    I(inds) = (kfreq-1)*2 + [1 1 2 2];
    J(inds) = (kfreq-1)*4 + [1 2 3 4];
    E(inds) = localMatrix_stator;
    
    %rotor surface
    localMatrix_stator = nu0*kfreq/r_si * [1 1;-1 -1];
    inds = (kfreq-1)*8 + (5:8);
    
    I(inds) = (kfreq-1)*2 + 2*Nfreqs + [1 1 2 2];
    J(inds) = (kfreq-1)*4 + [1 2 3 4];
    E(inds) = localMatrix_stator;
end
ag2H = sparse(I, J, E, 4*Nfreqs, 4*Nfreqs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling the restriction matrix Q: a matrix for selecting the air-gap
% boundary potentials from amongst all the potentials 
% [As; Ar] = Q*A;

Q = sparse(1:(N_ag_s + N_ag_r), bndData.agNodes_global, ones(1, N_ag_s + N_ag_r), N_ag_s + N_ag_r, Np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling the Neumann boundary matrix N: a matrix with the entries
% N_{ij} = Int_bnd -Phi_i*Ht_j
% with the rows 1...N_ag_s corresponding to the stator-ag boundary;
% and the remaining rows similarly for rotor

N = [r_si*MF(:, 1:N_ag_s)';
    -r_ro*MF(:, (N_ag_s + 1):(N_ag_s + N_ag_r))'];

end

