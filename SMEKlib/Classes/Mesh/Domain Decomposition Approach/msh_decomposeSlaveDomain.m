function msh = msh_decomposeSlaveDomain(msh, dims, els)
%msh_decomposeSlaveDomain Remove slave domain from mesh etc.
% 
% (c) 2018 Antti Lehikoinen / Aalto University

p_orig = msh.p;
t_orig = msh.t;

%removing domain from mesh:
msh = msh_removeElements(msh, cell2mat(els));
Np = size(msh.p, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing slave domain
%finding boundaries and defining slot polygon
[n_out, e_out] = findOuterNodes(t_orig(:, horzcat(els{1})));

rotA = pi/dims.Qs;
Xslot = [cos(-rotA) -sin(-rotA);sin(-rotA) cos(-rotA)]*p_orig(:,n_out);

rc = dims.rc;

%packing slot
figure(5); clf;
plot( Xslot(1,[1:end 1]), Xslot(2,[1:end 1]), 'b-'); axis equal;
[Xc, reff] = pack_slot_hexagonal(Xslot, dims.Nc_slot, rc);

%final adjustment
Xc(1,:) = Xc(1,:) - 0.0e-3;
msh.misc.Xc = Xc;

%msh.misc.Xc = Xc;
%msh.misc.Xslot = Xslot;

%taking the inter-layer insulation into account in a reaaaally ugly way
%Xc(1, (dims.Nc_slot/2 + 1):end) = Xc(1, (dims.Nc_slot/2 + 1):end) - 1e-3;
drawConductors(reff, Xc, 'EdgeColor', 'r'); axis equal;

% generating mesh for slot

%now with this fancy new pdetoolbox
slot_model = createpde(1);

%geometry description
slot_g_bnd = [2 size(Xslot,2) Xslot(1,:) Xslot(2,:)]'; %slot boundaries
slot_g_strands = [ones(1, dims.Nc_slot); Xc; rc*ones(dims.Nc_slot)];

slot_g = zeropadcat(slot_g_bnd, slot_g_strands(:,:));

%namespace
slot_ns = ['b' char('s'*ones(1,dims.Nc_slot))];

%set formula
slot_sf = 'b';

%decomposed geometry object
[slot_dl, slot_bt] = decsg(slot_g, slot_sf, slot_ns);

%adding geometry to the model
geometryFromEdges(slot_model, slot_dl);

%meshing
Hmax = rc/2;
if Hmax > 0
    generateMesh(slot_model, 'GeometricOrder', 'linear', 'Hmax', Hmax);
else
    generateMesh(slot_model, 'GeometricOrder', 'linear');
end
%
msh_slave = SimpleMesh(slot_model.Mesh.Nodes, slot_model.Mesh.Elements);

%getting conductor elements
[~, ~, t_temp] = meshToPet(slot_model.Mesh);
conductors_slave = get_elementsInDomain(slot_bt(:,2:end), t_temp(4,:));

%plotting
figure(6); clf; hold on; box on;
%pdeplot(slot_model); axis equal;
msh_triplot(msh_slave, [], 'b'); axis equal;
msh_triplot(msh_slave, horzcat(conductors_slave{:}), 'r');
msh_triplot(msh_slave, conductors_slave{1}, 'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating interpolation data from master to slave mesh
Np_slave = size(msh_slave.p, 2);
nd_slave = findOuterNodes(msh_slave.t);

msh.misc.msh_slave = msh_slave;
msh.misc.nd_slave = nd_slave;
msh.misc.conductors_slave = conductors_slave;

%interpolation from the boundary nodes of the removed mesh part, to the
%slave mesh
X_bnd = Xslot; 
N_Dir =size(X_bnd, 2);
e_bnd = [1:N_Dir;2:N_Dir 1];
m_c2s = generateInterpolants(X_bnd, msh_slave.p(:,nd_slave), e_bnd, 1:numel(nd_slave), @(X)( X ) );
P_c2s = sparse(m_c2s(:,1), m_c2s(:,2), m_c2s(:,3), numel(nd_slave), N_Dir);

%mapping master to slave nodes, plus periodic boundaries
Qs_sector = dims.Qs / msh.symmetrySectors;
P_data_master = cell(Qs_sector, 1);

figure(7); clf; hold on; msh_triplot(msh, [], 'k');
for kslot = 1:Qs_sector
    rotA = 2*pi/dims.Qs*(0.5 + kslot-1);
    rotM = [cos(rotA) -sin(rotA);sin(rotA) cos(rotA)];
    
    %finding master mesh nodes on the boundary
    X_bnd = rotM * Xslot;
    Idx = knnsearch(msh.p', X_bnd');
    
    %skipping interpolation since we have a 1-on-1 correspondence
    P_data_master{kslot} = sparse(1:N_Dir, Idx, ones(1, N_Dir), N_Dir, Np);
    %P_data_master{1, 2 + kslot} = [Idx'; (1:N_Dir) + Np_master + (kslot-1)*N_Dir; ones(1, N_Dir)];
    
    figure(7); plot(msh.p(1,Idx), msh.p(2,Idx), 'ro-')
    plot(X_bnd(1,:), X_bnd(2,:), 'k.')
end
msh.misc.P_D2s = P_c2s;
msh.misc.P_data_master = P_data_master;
msh.misc.P_m2D = cell2mat( P_data_master );

end