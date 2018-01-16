function [bndData, msh] = initialize_AGEboundaryData(msh, statorElements, rotorElements, arg_1)

if isa(arg_1, 'struct')
    %error('Not yet impemented!');
    dims = arg_1;
    agElements = [];
    Dp = 2*sum(msh.p.^2,1).^0.5;
    agNodes_all = [find( abs(Dp - dims.D_si) < 5e-4 ) ...
        find( abs(Dp - dims.D_ro) < 5e-4 )];
elseif size(arg_1, 1) == 1
    agElements = arg_1;
    t_ag = msh.t(:, agElements);
    agNodes_all = unique(t_ag(:))';
elseif size(arg_1, 1) == 3
    t_ag = arg_1;
    agNodes_all = unique(t_ag(:))';
    agElements = [];
end

%identifying nodes
statorNodes_all = unique( msh.t(:,statorElements(:)) )';
rotorNodes_all = unique( msh.t(:, rotorElements(:)) )';

%ag-nodes on stator and rotor surface
agNodes_rotor = intersect(rotorNodes_all, agNodes_all);
agNodes_stator = intersect(statorNodes_all, agNodes_all);
agNodes_middle = setdiff(agNodes_all, union(agNodes_rotor, agNodes_stator));

Np_old = size(msh.p, 2);

%removing middle nodes
[~,IA] = setdiff(1:Np_old, agNodes_middle);
newInds = zeros(1, Np_old); newInds(IA) = 1:numel(IA);

msh.p = msh.p(:, IA);
msh.t = msh.t(:, setdiff(1:size(msh.t,2), agElements));
msh.t = reshape(newInds( msh.t(:)), 3, []);

agNodes_rotor = newInds( agNodes_rotor );
agNodes_stator = newInds( agNodes_stator );

N_ag_r = numel(agNodes_rotor);
N_ag_s = numel(agNodes_stator);

%sorting rotor and stator nodes
agAngles_rotor = atan2( msh.p(2, agNodes_rotor), msh.p(1, agNodes_rotor) );
agAngles_rotor( agAngles_rotor < 0 ) = agAngles_rotor( agAngles_rotor < 0 ) + 2*pi;
[agAngles_rotor, agOrder_rotor] = sort( agAngles_rotor );

agAngles_stator = atan2( msh.p(2, agNodes_stator), msh.p(1, agNodes_stator) );
agAngles_stator( agAngles_stator < 0 ) = agAngles_stator( agAngles_stator < 0 ) + 2*pi;
[agAngles_stator, agOrder_stator] = sort( agAngles_stator );

agAngles_global = [agAngles_stator agAngles_rotor];
agNodes_global = [agNodes_stator(agOrder_stator) agNodes_rotor(agOrder_rotor)];

%symmetry sectors --> generating virtual data
if isfield(msh, 'symmetrySectors')
    symm = msh.symmetrySectors;
    
    %virtual air-gap nodes and the symmetry sectors they belong to
    ps_temp = [msh.p(:, agNodes_stator) zeros(2, (symm-1)*N_ag_s)]; sec_s = zeros(1, symm*N_ag_s);
    pr_temp = [msh.p(:, agNodes_rotor) zeros(2, (symm-1)*N_ag_r)]; sec_r = zeros(1, symm*N_ag_r);
    
    for k_sector = 1:(symm-1)
        %replicating nodes
        rA = k_sector * 2*pi/symm; rotM = [cos(rA) -sin(rA); sin(rA) cos(rA)];
        ps_temp(:, k_sector*N_ag_s + (1:N_ag_s)) = rotM * msh.p(:, agNodes_stator);
        pr_temp(:, k_sector*N_ag_r + (1:N_ag_r)) = rotM * msh.p(:, agNodes_rotor);
        
        sec_s(1, k_sector*N_ag_s + (1:N_ag_s)) = k_sector;
        sec_r(1, k_sector*N_ag_r + (1:N_ag_r)) = k_sector;
    end
    
    %identities of the virtual nodes
    Is = repmat(1:numel(agNodes_stator), 1, symm);
    Ir = repmat(1:numel(agNodes_rotor), 1, symm);
    
    %removing duplicate nodes on the boundaries
    ind_s = setdiff( 1:(symm*N_ag_s), [(N_ag_s+1):N_ag_s:(symm*N_ag_s) (symm*N_ag_s)] );
    ind_r = setdiff( 1:(symm*N_ag_r), [(N_ag_r+1):N_ag_r:(symm*N_ag_r) (symm*N_ag_r)] );
    
    N_ag_s_true = N_ag_s;
    N_ag_r_true = N_ag_r;
    
    N_ag_s = numel(ind_s);
    N_ag_r = numel(ind_r);
    
    p_ag_virt = [ps_temp(:, ind_s) pr_temp(:, ind_r)];
    
    el_table = [1:size(p_ag_virt,2); Is(ind_s) N_ag_s_true+Ir(ind_r); sec_s(ind_s) sec_r(ind_r)];
    
    %re-sorting rotor and stator nodes
    agAngles_rotor = atan2( p_ag_virt(2, (N_ag_s+1):end), p_ag_virt(1, (N_ag_s+1):end) );
    agAngles_rotor( agAngles_rotor < 0 ) = agAngles_rotor( agAngles_rotor < 0 ) + 2*pi;
    agAngles_rotor = sort( agAngles_rotor );

    agAngles_stator = atan2( p_ag_virt(2, 1:N_ag_s), p_ag_virt(1, 1:N_ag_s) );
    agAngles_stator( agAngles_stator < 0 ) = agAngles_stator( agAngles_stator < 0 ) + 2*pi;
    agAngles_stator = sort( agAngles_stator );

    agAngles_global = [agAngles_stator agAngles_rotor];
end

bndData = struct('N_ag_s', N_ag_s, 'N_ag_r', N_ag_r, ...
    'agNodes_global', agNodes_global, 'agAngles_global', agAngles_global);

if isfield(msh, 'symmetrySectors')
    bndData.p_ag_virt = p_ag_virt;
    bndData.el_table = el_table;
    bndData.N_ag_s_true = N_ag_s_true;
    bndData.N_ag_r_true = N_ag_r_true;
end

end