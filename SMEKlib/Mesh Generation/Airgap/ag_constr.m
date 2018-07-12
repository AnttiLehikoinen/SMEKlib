function this = ag_constr(this, msh, varargin)
%

%FIXME: add number of midnodes thingy

% generating air-gap triangulation; for two layers by default
this.misc = struct();

if numel(varargin) && isa(varargin{1}, 'double') && size(varargin{1}, 1) > 1
    % triangulation given as input
    [n_ag_s, n_ag_r, N_ag_s, N_ag_r, n_mid, t_in, t_out] = givenTriangulation(msh, varargin{1});
    N_mid = numel(n_mid);
    mid_given = true;
else
    %generating triangulatio
    n_ag_s = msh.namedNodes.get('n_ag_s');
    n_ag_r = msh.namedNodes.get('n_ag_r');
    N_ag_s = numel(n_ag_s);
    N_ag_r = numel(n_ag_r);

    %sorting
    %[~,I] = sort( atan2(msh.p(2,n_ag_s), msh.p(1,n_ag_s)) ); n_ag_s = n_ag_s(I);
    %[~,I] = sort( atan2(msh.p(2,n_ag_r), msh.p(1,n_ag_r)) ); n_ag_r = n_ag_r(I);

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
        %N_mid = 2*(N_ag_s-1) + 1;
        %N_mid = ceil( 0.9*N_ag_r )

        %angles = linspace(0, 2*pi/msh.symmetrySectors, N_mid);
        angle_start = (atan2(msh.p(2,n_ag_s(1)), msh.p(1,n_ag_s(1))) + ...
            atan2(msh.p(2,n_ag_r(1)), msh.p(1,n_ag_r(1))))/2;
        angles = linspace(0, 2*pi/msh.symmetrySectors, N_mid) + angle_start;

        rin = mean( sum(msh.p(:,n_ag_r).^2,1).^0.5 );
        rout = mean( sum(msh.p(:,n_ag_s).^2,1).^0.5 );

        p_mid = (rin+rout)/2*[cos(angles); sin(angles)];
        n_mid = size(msh.p, 2) + (1:N_mid);

        %adding mid-nodes to the mesh
        msh.p = [msh.p p_mid];
    end

    %layer-triangulations
    t_in = singleLayerAGtriangulation_2(msh, n_mid, n_ag_r);
    t_out = singleLayerAGtriangulation_2(msh, n_ag_s, n_mid);
end

agNodes_global = [n_ag_s n_ag_r n_mid];

%{
figure(12); clf; hold on; axis equal;
p = msh.p;
triplot(t_in(:,:)', p(1,:), p(2,:), 'r');
%t_in(:,1:10)
triplot(t_out', p(1,:), p(2,:), 'b');

plot(p(1, n_mid), p(2,n_mid), 'ro-');
plot(p(1, n_ag_s), p(2, n_ag_s), 'bo-');
plot(p(1, n_ag_r), p(2, n_ag_r), 'go-');
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
if symm > 1
    Ntemp = size(p_ag_virt, 2);
    inds2keep = setdiff(1:Ntemp, [((N_ag_s+N_ag_r + N_mid +1):N_mid:Ntemp) Ntemp]);
    p_ag_virt = p_ag_virt(:, inds2keep);
    virt_sectors = virt_sectors(inds2keep);
    virt_identities = virt_identities(inds2keep);
end

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

%setting constant part of AG matrix
this.setConstantAGmatrix(size(msh.p,2));

this.msh_ag.t = this.t_moving;

end



function [n_ag_s, n_ag_r, N_ag_s, N_ag_r, n_mid, t_in, t_out] = givenTriangulation(msh, tag)

rotel = msh.rotel;
statel = setdiff(1:size(msh.t,2), rotel);

n_ag_s = toRow(unique(intersect( msh.t(:,statel), tag)));
n_ag_r = toRow(unique(intersect( msh.t(:,rotel), tag)));
n_mid = toRow(setdiff(tag, [n_ag_s n_ag_r]));

%sorting
n_ag_s = nsort(msh, n_ag_s);
n_ag_r = nsort(msh, n_ag_r);
n_mid = nsort(msh, n_mid);

msh.namedNodes.add('n_ag_s', n_ag_s);
msh.namedNodes.add('n_ag_r', n_ag_r);
msh.namedNodes.add('n_mid', n_mid);

N_ag_s = numel(n_ag_s);
N_ag_r = numel(n_ag_r);

inds_out = find( sum( ismember(tag, n_ag_s), 1) );
t_out = tag(:,inds_out);
t_in = tag(:, setdiff(1:size(tag,2), inds_out));

end

function n = nsort(msh, n)

atemp = atan2(msh.p(2,n), msh.p(1,n)); atemp( atemp<0 ) = atemp(atemp<0) + 2*pi;
[~, I] =sort(atemp);
n = n(I);

end



