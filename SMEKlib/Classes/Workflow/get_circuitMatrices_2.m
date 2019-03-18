function [Sc, Mc, Cstar2, Ustar] = get_circuitMatrices_2(sim, varargin)
%get_circuitMatrices returns all circuit matrices for simulation.
%
% Now without elimination of rotor cage currents.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

Cstar = [1 0; 0 1;-1 -1];
Ustar = eye(3); Cstar2 = [1 0; 0 1;-1 -1];
if isfield(sim.dims, 'N_inParallel')
    n = sim.dims.N_inParallel;
    Mtemp = zeros(n, n-1);
    for k = 1:(n-1)
        Mtemp([k k+1], k) = [1;-1];
    end
    
    Cstar2 = [kron(Cstar, eye(n)) ...
        [zeros(2*n, n-1); Mtemp]];
    Ustar = [kron(eye(2), ones(n,1)); zeros(n-1, 2)];
end

Np = sim.Np;
if sim.dims.type_statorWinding == defs.stranded || sim.dims.type_statorWinding == defs.decomposed
    Ni_s = size(sim.matrices.Ls, 2); Nu_s = 0;
    S_AU_s = sparse(Np, 0);
    S_AI_s = -sim.matrices.Cs*sim.matrices.Ls;
    
    M_UA_s = sparse(0, Np);
    
    S_UU_s = [];
    M_UU_s = [];
    
    S_UI_s = sparse(Nu_s, Ni_s);
    
    M_IA_s = sim.dims.leff*sim.matrices.Ls'*sim.matrices.Cs';
    S_IU_s = sparse(Ni_s, 0);
    
    S_II_s = real(sim.matrices.Zew_s);
    M_II_s = imag(sim.matrices.Zew_s);
end
if sim.dims.type_statorWinding == defs.decomposed
   S_AI_s = S_AI_s*0;
   M_IA_s = M_IA_s*0;
end

%checking if star connection
if sim.dims.connection_stator == defs.star && sim.dims.type_statorWinding ~= defs.decomposed
   S_AI_s = S_AI_s*Cstar;
   S_UI_s = S_UI_s*Cstar;
   M_IA_s = Cstar'*M_IA_s;
   S_IU_s = Cstar'*S_IU_s;
   S_II_s = Cstar'*S_II_s*Cstar;
   M_II_s = Cstar'*M_II_s*Cstar;
   Ni_s = 2;
elseif sim.dims.connection_stator == defs.star && sim.dims.type_statorWinding == defs.decomposed
   S_AI_s = S_AI_s*Cstar2;
   S_UI_s = S_UI_s*Cstar2;
   M_IA_s = Cstar2'*M_IA_s;
   S_IU_s = Cstar2'*S_IU_s;
   S_II_s = Cstar2'*S_II_s*Cstar2;
   M_II_s = Cstar2'*M_II_s*Cstar2;
   Ni_s = size(Cstar2, 2);
end
    

if isfield(sim.dims, 'supply_type') && sim.dims.supply_type == defs.current_supply
    S_AI_s = sparse(Np, 0);
    S_AU_s = sparse(Np, 0);
    S_UI_s = [];
    M_IA_s = sparse(0, Np);
    S_IU_s = [];
    S_II_s = [];
    M_II_s = [];
    Ni_s = 0; Nu_s = 0;
elseif isfield(sim.dims, 'supply_type') && sim.dims.supply_type == defs.current_supply_dynamic
    Nphases = 3;
    
    Nu_temp = size(S_UI_s,1);
    S_AI_s = [S_AI_s sparse(Np, Nphases)];
    S_UI_s = [S_UI_s sparse(Nu_temp, Nphases)];
    M_IA_s = [M_IA_s; sparse(Nphases, Np)];
    S_IU_s = [S_IU_s; sparse(Nphases, Nu_temp)];
    
    S_II_s = [S_II_s -kron( eye(Nphases), ones(Ni_s/Nphases, 1) );
        kron( eye(Nphases), ones(Ni_s/Nphases, 1) )' sparse(Nphases,Nphases)];
    %full(S_II_s)
    M_II_s = [M_II_s sparse(Ni_s, Nphases);
        sparse(Nphases, Ni_s+Nphases)];
    
    %{
    S_II_s = [S_II_s -speye(Ni_s, Ni_s);
        speye(Ni_s,Ni_s) sparse(Ni_s,Ni_s)];
    M_II_s = [M_II_s sparse(Ni_s, Ni_s);
        sparse(Ni_s, 2*Ni_s)];
    %}
    Ni_s = Ni_s + Nphases;
end
    
   
sim.results.Ni_s = Ni_s; sim.results.Nu_s = Nu_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim.dims.type_rotorWinding == defs.cage || sim.dims.type_rotorWinding == defs.user_defined
    Nu_r = size(sim.matrices.Lr, 1); 
    Ni_r = size(sim.matrices.Lr, 2);
    
    S_AU_r = -1/sim.dims.leff*sim.matrices.Cr;
    S_AI_r = sparse(Np, Ni_r);
    
    M_UA_r = sim.matrices.DRr*transpose(sim.matrices.Cr);
    
    S_UU_r = -speye(Nu_r);
    M_UU_r = sparse(Nu_r, Nu_r);
    
    S_UI_r = sim.matrices.DRr*sim.matrices.Lr;
    
    M_IA_r = sparse(Ni_r, Np);
    S_IU_r = sim.matrices.Lr';
    
    S_II_r = real(sim.matrices.Zew_r);
    M_II_r = imag(sim.matrices.Zew_r);
else
    %no rotor circuit
    S_AU_r = []; S_AI_r = [];
    
    M_UA_r = []; S_UU_r = []; M_UU_r = [];
    
    S_UI_r = []; M_IA_r = []; S_IU_r = [];
    S_II_r = []; M_II_r = [];
    Nu_r = 0; Ni_r = 0;
end     

sim.results.Nu_r = Nu_r; sim.results.Ni_r = Ni_r;

if numel(varargin)
    rc = varargin{1}; %slip
else
    rc = 1;
end

Nu = Nu_s + Nu_r;
Ni = Ni_s + Ni_r;

Sc = [sparse(Np, Np) S_AU_s S_AU_r S_AI_s S_AI_r;
    sparse(Nu, Np) blkdiag(S_UU_s, S_UU_r) blkdiag(S_UI_s, S_UI_r);
    sparse(Ni, Np) blkdiag(S_IU_s, S_IU_r) blkdiag(S_II_s, S_II_r)];

Mc = [sim.matrices.Ms+rc*sim.matrices.Mr sparse(Np, Nu+Ni);
    [M_UA_s; rc*M_UA_r] blkdiag(M_UU_s, M_UU_r) sparse(Nu, Ni);
    [M_IA_s; rc*M_IA_r] sparse(Ni, Nu) blkdiag(M_II_s, rc*M_II_r)];

%{
mats = struct('S_AU', [S_AU_s S_AU_r], 'S_AI', [S_AI_s S_AI_r], ...
    'M_UA', [M_UA_s; M_UA_r], ...
    'S_UU', blkdiag(S_UU_s, S_UU_r), 'M_UU', blkdiag(M_UU_s, M_UU_r), ...
    'S_UI', [S_UI_s; S_UI_r], ...
    'M_IA', [M_IA_s; M_IA_s], 'S_IU',
%}