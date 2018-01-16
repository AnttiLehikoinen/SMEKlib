function [bndData, msh] = initialize_slidingBoundaryData(msh, statorElements, rotorElements, varargin)

if numel(varargin) == 0
    error('Not yet impemented!');
else
    agElements = varargin{1};
    t_ag = msh.t(:, agElements);
end

%identifying nodes

statorNodes_all = unique( msh.t(:,statorElements(:)) )';
rotorNodes_all = unique( msh.t(:, rotorElements(:)) )';
agNodes_all = unique(t_ag(:))';

%ag-nodes on stator and rotor surface
agNodes_rotor = intersect(rotorNodes_all, agNodes_all);
agNodes_stator = intersect(statorNodes_all, agNodes_all);

agNodes_middle = setdiff(agNodes_all, union(agNodes_rotor, agNodes_stator));

if ~any(agNodes_middle)
    error('Single-band functionality not yet implemented!');
end

Np_old = size(msh.p, 2); N_ag_m = numel(agNodes_middle);
msh.p = [msh.p msh.p(:, agNodes_middle)];
temp_inds = 1:Np_old; temp_inds( agNodes_middle ) = (Np_old + 1):(Np_old + N_ag_m);

%air-gap elements based on the rotor
el_ag_rotorside = agElements(find( sum(ismember(t_ag, agNodes_rotor),1) > 0));

msh.t(:, el_ag_rotorside) = reshape( temp_inds( msh.t(:,el_ag_rotorside) ), 3, []);

%changing the moving surface
agNodes_rotor = (Np_old + 1):(Np_old + N_ag_m);
agNodes_stator = agNodes_middle;
N_ag_r = N_ag_m;
N_ag_s = N_ag_m;

%sorting rotor and stator nodes
agAngles_rotor = atan2( msh.p(2, agNodes_rotor), msh.p(1, agNodes_rotor) );
[agAngles_rotor, agOrder_rotor] = sort( agAngles_rotor );

agAngles_stator = atan2( msh.p(2, agNodes_stator), msh.p(1, agNodes_stator) );
[agAngles_stator, agOrder_stator] = sort( agAngles_stator );

agAngles_global = [agAngles_stator agAngles_rotor];
agNodes_global = [agNodes_stator(agOrder_stator) agNodes_rotor(agOrder_rotor)];

bndData = struct('N_ag_s', N_ag_s, 'N_ag_r', N_ag_r, ...
    'agNodes_global', agNodes_global, 'agAngles_global', agAngles_global);

end