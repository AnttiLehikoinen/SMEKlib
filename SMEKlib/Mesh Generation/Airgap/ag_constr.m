function this = ag_constr(this, msh, varargin)

% generating air-gap triangulation; for two layers by default
this.misc = struct();

n_ag_s = msh.namedNodes.get('n_ag_s');
n_ag_r = msh.namedNodes.get('n_ag_r');
N_ag_s = numel(n_ag_s);
N_ag_r = numel(n_ag_r);

%sorting
[~,I] = sort( atan2(msh.p(2,n_ag_s), msh.p(1,n_ag_s)) ); n_ag_s = n_ag_s(I);
[~,I] = sort( atan2(msh.p(2,n_ag_r), msh.p(1,n_ag_r)) ); n_ag_r = n_ag_r(I);

n_mid = msh.namedNodes.get('n_mid');
mid_given = false;
if any(n_mid)
    N_mid = numel(n_mid);
    [~,I] = sort( atan2(msh.p(2,n_mid), msh.p(1,n_mid)) ); n_mid = n_mid(I);
    
    p_mid = msh.p(n_mid);
    mid_given = true;
else
    %no middle layer nodes given --> generating middle layer
    N_mid = ceil( 0.5*(N_ag_s + N_ag_r) ) - 1;
    %N_mid = ceil(0.5*N_ag_s) - 1;
    N_mid = 2*(N_ag_s-1) + 1;

    angles = linspace(0, 2*pi/msh.symmetrySectors, N_mid);
    rin = mean( sum(msh.p(:,n_ag_r).^2,1).^0.5 );
    rout = mean( sum(msh.p(:,n_ag_s).^2,1).^0.5 );

    p_mid = (rin+rout)/2*[cos(angles); sin(angles)];
    n_mid = size(msh.p, 2) + (1:N_mid);

    %adding mid-nodes to the mesh
    msh.p = [msh.p p_mid];
end

agNodes_global = [n_ag_s n_ag_r n_mid];

%layer-triangulations
t_in = singleLayerAGtriangulation_2(msh, n_mid, n_ag_r);
t_out = singleLayerAGtriangulation_2(msh, n_ag_s, n_mid);

%{
figure(12); clf; hold on;
p = msh.p;
plot(p(1, n_mid), p(2,n_mid), 'ro-');
triplot(t_in(:,1:10)', p(1,:), p(2,:), 'r');
t_in(:,1:10)
triplot(t_out', p(1,:), p(2,:), 'b');
%}

%switching from global to local indexing
[~, inds_local] = ismember(t_in(:), agNodes_global);
t_in = reshape(inds_local(:), size(t_in,1), []);
[~, inds_local] = ismember(t_out(:), agNodes_global);
t_out = reshape(inds_local(:), size(t_out,1), []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating stuff necessary for movement

% generating virtual nodes
symm = msh.symmetrySectors;

%adding virtual rotor surface nodes
p_ag_virt = [msh.p(:,agNodes_global) zeros(2, (symm-1)*N_mid)];
virt_sectors = zeros(1, size(p_ag_virt, 2));
virt_identities = [n_ag_s n_ag_r repmat(n_mid, 1, symm)];
for k_sector = 1:(symm-1)
    rA = k_sector * 2*pi/symm;
    
    inds = (N_ag_s + N_ag_r + N_mid) + (1:N_mid) + (k_sector-1)*N_mid;
    
    p_ag_virt(:, inds) = ...
        [cos(rA) -sin(rA);sin(rA) cos(rA)] * msh.p(:, n_mid);

    virt_sectors(inds) = k_sector;
end

%getting rid of redundant dublicate virtual nodes
Ntemp = size(p_ag_virt, 2);
inds2keep = setdiff(1:Ntemp, [((N_ag_s+N_ag_r + N_mid +1):N_mid:Ntemp) Ntemp]);
p_ag_virt = p_ag_virt(:, inds2keep);
virt_sectors = virt_sectors(inds2keep);
virt_identities = virt_identities(inds2keep);

% indices of moving interface node within the moving-band triangulation
inds_r = find( ismember(agNodes_global(t_out(:)), n_mid ) );
sortedNodes_moving = (N_ag_s + N_ag_r + 1):size(p_ag_virt,2);
[~, originalPositions_rotor] = ismember( t_out(inds_r), sortedNodes_moving );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving data
this.agNodes_global = agNodes_global;
this.t_const = t_in;
this.p_virt = p_ag_virt;

%this.n_moving = [n_ag_r n_mid];
this.t_moving = t_out;
this.n_moving = (N_ag_s + 1):size(p_ag_virt,2);
this.n_bnd = (N_ag_s + N_ag_r + 1):size(p_ag_virt,2);

this.inds_r = inds_r;
this.original_positions = originalPositions_rotor;

this.shiftTol = 2*pi / numel(this.n_bnd);
%this.n_moving_local = sortedNodes_moving;

this.msh_ag = MachineMesh(this.p_virt, this.t_moving);

this.el_table = [1:size(p_ag_virt, 2);
            virt_identities;
            msh.periodicityCoeff.^virt_sectors];
        
%adding periodic nodes
if ~mid_given
    msh.namedNodes.add('Periodic_master', n_mid(1));
    msh.namedNodes.add('Periodic_slave', n_mid(end));
end
        
%generating constant-part of ag matrix
Np = size(msh.p,2);
this.msh_ag.t = this.t_const;
Sag_c = ...
    MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), 1/(pi*4e-7), [], this.msh_ag);
%moving to global indexing and taking care of symmetry sectors
inds = 1:Sag_c.Nvals;
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.I(inds));
Sag_c.I(inds) = this.el_table(2, Sag_c.I(inds));
Sag_c.E(inds) = Sag_c.E(inds) .* this.el_table(3, Sag_c.J(inds));
Sag_c.J(inds) = this.el_table(2, Sag_c.J(inds));
this.S_const = Sag_c.finalize(Np,Np);

this.msh_ag.t = this.t_moving;

end

