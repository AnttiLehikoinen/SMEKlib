function bandData = initializeBandData(msh, statorElements, rotorElements, varargin)
%initializeBandData initializes moving band data.
%
% bandData = initializeBandData(msh, statorElements, rotorElements, dims, t_ag)
%   computes the bandData struct from the given mesh and the air-gap
%   triangulation t_ag (using global indexing to nodes)
% bandData = initializeBandData(msh, statorElements, rotorElements, dims)
%   tries to generate the air-gap triangulation based on the machine
%   dimensions
%
% The bandData struct contains the following fields related to the air-gap
% triangulation:
%   - Ne_ag, N_ag_s, N_ag_r: number of air-gap elements, and static and
%   moving nodes respectively.
%   - t_ag: air-gap triangulation, with local indexing 1:(N_ag_s+N_ag_r)
%   - agNodes_global: indices for switching from local to global indexing:
%   inds_global = agNodes_global(inds_local)
%   - agAngles_all: angular coordinates of the air-gap nodes
%   - shiftTol: period after which t_ag changes with rotation
%   - sortedNodes_rotor: moving nodes in angle-sorted order, local indexing
%   - inds_r: linear indices of moving nodes inside t_ag(:)
%   - originalPositions_rotor: original indices of t_ag(inds_r) inside
%   sortedNodes_rotor
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if numel(varargin) == 1
    t_ag = varargin{1};
    
    rotorNodes_all = unique( msh.t(:, rotorElements(:)) )';
    agNodes_all = unique(t_ag(:))';
    
    agNodes_rotor = intersect(rotorNodes_all, agNodes_all);
    agNodes_stator = setdiff(agNodes_all, rotorNodes_all);
   
    %un-used stuff for two-band movement
    %statorNodes_all = unique( msh.t(:,statorElements(:)) )';
    %agNodes_middle = setdiff(agNodes_all, union(intersect(statorNodes_all, agNodes_all), agNodes_rotor));
    %r_all = sum(msh.p(:, agNodes_all).^2,1).^0.5;
    %middleNode_shiftRatio = (r_all - dims.D_ro/2) / ((dims.D_si - dims.D_ro)/2);
elseif numel(varargin) == 0
    error('Not yet implemented');
else
    error('Incorrect number of input arguments!')
end

Ne_ag = size(t_ag, 2);
N_ag_s = numel( agNodes_stator );
N_ag_r = numel( agNodes_rotor );

%sorting rotor nodes
agAngles_rotor = atan2( msh.p(2, agNodes_rotor), msh.p(1, agNodes_rotor) );
agAngles_rotor( agAngles_rotor<0 ) = agAngles_rotor(agAngles_rotor<0) + 2*pi;
[~, agOrder_rotor] = sort( agAngles_rotor );
agNodes_rotor = agNodes_rotor(agOrder_rotor);

agNodes_global = [agNodes_stator agNodes_rotor];

%switching from global to local indexing
[~, inds_local] = ismember(t_ag(:), agNodes_global);
t_ag = reshape(inds_local(:), 3, []);



%sectorAngle = max(agAngles_rotor) - min(agAngles_rotor);

if isfield(msh, 'symmetrySectors')
    % periodicity conditions --> adding virtual nodes on the ag-rotor
    % surface to encompass the entire cross-section  
    symm = msh.symmetrySectors;
    
    %TODO: sort nodes first!
    
    %adding virtual rotor surface nodes
    p_ag_virt = [msh.p(:,agNodes_global) zeros(2, (symm-1)*N_ag_r)];
    virt_sectors = zeros(1, size(p_ag_virt, 2));
    virt_identities = [agNodes_stator repmat(agNodes_rotor, 1, symm)];
    for k_sector = 1:(symm-1)
        rA = k_sector * 2*pi/symm;
        p_ag_virt(:, (N_ag_s+N_ag_r)+(1:N_ag_r) + (k_sector-1)*N_ag_r) = ...
            [cos(rA) -sin(rA);sin(rA) cos(rA)] * msh.p(:, agNodes_rotor);
        
        virt_sectors((N_ag_s+N_ag_r)+(1:N_ag_r) + (k_sector-1)*N_ag_r) = k_sector;
    end
    
    %getting rid of redundant dublicate rotor nodes
    Ntemp = size(p_ag_virt, 2);
    p_ag_virt = p_ag_virt(:, ...
        setdiff(1:Ntemp, [((N_ag_s+N_ag_r+1):N_ag_r:Ntemp) Ntemp]));
    virt_sectors = virt_sectors(...
        setdiff(1:Ntemp, [((N_ag_s+N_ag_r+1):N_ag_r:Ntemp) Ntemp]));
    virt_identities = virt_identities(...
        setdiff(1:Ntemp, [((N_ag_s+N_ag_r+1):N_ag_r:Ntemp) Ntemp]));
    
    N_ag_r = size(p_ag_virt, 2) - N_ag_s;
   
    
    %sorting rotor and stator nodes
    agAngles_rotor = atan2( p_ag_virt(2,(N_ag_s+1):end), p_ag_virt(1,(N_ag_s+1):end) );
    agAngles_rotor( agAngles_rotor<0 ) = agAngles_rotor(agAngles_rotor<0) + 2*pi;
    [~, agOrder_rotor] = sort( agAngles_rotor );

    %sectorAngle = max(agAngles_rotor) - min(agAngles_rotor)
    sectorAngle = 2*pi;
else
    %sorting rotor and stator nodes
    %agAngles_rotor = atan2( msh.p(2, agNodes_rotor), msh.p(1, agNodes_rotor) );
    %agAngles_rotor( agAngles_rotor<0 ) = agAngles_rotor(agAngles_rotor<0) + 2*pi;
    %[~, agOrder_rotor] = sort( agAngles_rotor );

    %sectorAngle = max(agAngles_rotor) - min(agAngles_rotor)
    sectorAngle = 2*pi;
end

agAngles_stator = atan2( msh.p(2, agNodes_stator), msh.p(1, agNodes_stator) );
agAngles_stator( agAngles_stator<0 ) = agAngles_stator(agAngles_stator<0) + 2*pi;
%[~, agOrder_stator] = sort( agAngles_stator );

agAngles_all = [agAngles_stator agAngles_rotor];

%linear indices to rotor-surface nodes in t_ag
%inds_r = find( ismember(t_ag(:), agNodes_rotor ) );
inds_r = find( ismember(agNodes_global(t_ag(:)), agNodes_rotor ) );

%rotor nodes in sorted order; and the indices of the original air-gap
%triangulation in that list
temp = (N_ag_s + 1):(N_ag_s + N_ag_r); sortedNodes_rotor = temp(agOrder_rotor);
%sortedNodes_rotor = agNodes_rotor(agOrder_rotor);

[~, originalPositions_rotor] = ismember( t_ag(inds_r), sortedNodes_rotor );

%{
el_rotorBased = find( sum(ismember(t_ag, agNodes_rotor),1) == 2);
%finding rotor-based elements (= long edge on the rotor surface)
el_rotorBased = find( sum(ismember(t_ag, agNodes_rotor),1) == 2);

%sorting rotor elements
tempAngles = atan2( 0.5*(msh.p(2,t_ag(1,el_rotorBased))+msh.p(2,t_ag(3,el_rotorBased))), 0.5*(msh.p(1,t_ag(1,el_rotorBased))+msh.p(1,t_ag(3,el_rotorBased))) );
[~,rotorElementOrder] = sort( tempAngles );
el_rotorBased = el_rotorBased( rotorElementOrder );

%re-sorting?
elementOrder = 1:Ne_ag; Ne_ag_r = numel(el_rotorBased);
elementOrder(el_rotorBased) = el_rotorBased( mod((1:Ne_ag_r) - 1 - 0, Ne_ag_r) + 1 );
t_ag = t_ag(:, elementOrder); t_ag = sort(t_ag, 1);
%}


shiftTol = sectorAngle / N_ag_r;


bandData = struct('agNodes_global', agNodes_global, 'agAngles_all', agAngles_all, ...
    'Ne_ag', Ne_ag, 'N_ag_s', N_ag_s, 'N_ag_r', N_ag_r, ...
    't_ag', t_ag, 'shiftTol', shiftTol, ...
    'sortedNodes_rotor', sortedNodes_rotor, 'originalPositions_rotor', originalPositions_rotor, 'inds_r', inds_r);    

if isfield(msh, 'symmetrySectors')
    bandData.p_ag_virt = p_ag_virt;
    bandData.virt_sectors = virt_sectors;
    
    %temp = zeros(1, size(msh.p,2));
    %temp(agNodes_stator) = agNodes_stator;
    %temp(agNodes_rotor) = agNodes_rotor_orig;
    
    bandData.virt_identities = (virt_identities);
    
    bandData.el_table = [1:size(p_ag_virt, 2);
            bandData.virt_identities;
            msh.periodicityCoeff.^virt_sectors];
    
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting several non-necessary demonstrations here

plot(msh.p(1, agNodes_stator), msh.p(2, agNodes_stator), 'bo');
plot(msh.p(1, agNodes_rotor), msh.p(2, agNodes_rotor), 'ro');

%triplot(t_ag(:,el_rotorBased)', msh.p(1,:), msh.p(2,:), 'm')

%plotting a demonstration
figure(5); clf; hold on;
p_ag = msh.p; 
%plot( p_ag(1, t_ag(1,el_rotorBased)), p_ag(2, t_ag(1,el_rotorBased)), 'ro-')

plotmesh = msh; plotmesh.bandData = bandData;


Np = 20;
angles = linspace(0, 2.5, Np)/180*pi;
[~,n_plot] = min( abs(agAngles_rotor ) );
for kp = 1:Np
    rotorAngle = angles(kp);
    
    [t_ag_new, ~, p_ag] = updateRotorPosition(rotorAngle, plotmesh);
    t_ag_new = reshape( agNodes_global(t_ag_new(:)), 3, []);
    
    %{
    p_ag(:,rotorNodes_all) = [cos(rotorAngle) -sin(rotorAngle); sin(rotorAngle) cos(rotorAngle)] * msh.p(:,rotorNodes_all);
    %t_ag = updateRotorPosition(angles(kp), msh, bandData);

    nodeShift = -floor( (rotorAngle - 1*shiftTol/2) / shiftTol ) - 1;
    %newOrder = bandData.agOrder_rotor( mod((1:bandData.N_ag_r)-1 + nodeShift, bandData.N_ag_r) + 1 );
    newPositions = mod(originalPositions_rotor - 1 + nodeShift, N_ag_r ) + 1;
    
    %newPositions = agNodes_rotor( agOrder_rotor(newPositions) );
    %newPositions = sortedNodes_rotor(newPositions);
    
    t_ag_new = t_ag;
    %t_ag(inds_r) = agNodes_rotor(inds_agOrder_r(newPositions));
    t_ag_new(inds_r) = sortedNodes_rotor(newPositions);
    t_ag_new = reshape( agNodes_global(t_ag_new(:)), 3, []);
    %}
    
    figure(5); clf; hold on;
    triplot(t_ag_new', p_ag(1,:), p_ag(2,:), 'm');
    plot(p_ag(1, sortedNodes_rotor), p_ag(2,sortedNodes_rotor), 'ko-');
    plot(p_ag(1, agNodes_rotor(n_plot)), p_ag(2,agNodes_rotor(n_plot)), 'mx');
    axis([0.07 0.101 -0.02 0.02]);
    
    pause(0.2);
end


end