function S_ag = get_MovingBandMatrix_2(rotorAngle, varargin)
%get_MovingBandMatrix_2 moving band matrix with remeshing.
% 
% Like get_MovingBandMatrix, but now the air-gap is explicitly re-meshed.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

mu0 = pi*4e-7;
if (numel(varargin) == 1) || ( (numel(varargin)==2) && isa(varargin{2}, 'double'))
    msh = varargin{1};
    
    if ( (numel(varargin)==2) && isa(varargin{2}, 'double'))
        N = varargin{2};
    else
        N = size(msh.p,2);
    end
    
    if ~isfield(msh, 'symmetrySectors')
        error('Not implemented yet');
    else
        %number of nodes to skip on rotor surface
        nodeShift = -floor( (rotorAngle - 1*msh.bandData.shiftTol/2) / msh.bandData.shiftTol ) - 1;
        
        symm = msh.symmetrySectors;
        N_ag_r_b = (msh.bandData.N_ag_r + symm) / symm; %number of rotor nodes in symmetry sector
        N_ag_s = msh.bandData.N_ag_s;
        
        el_r = unique( floor( (msh.bandData.inds_r-1)/3 ) + 1); %rotor-bordering elements
        n_s = msh.bandData.t_ag(:, el_r); n_s = unique( n_s( n_s <= N_ag_s ) )'; %stator nodes connected to these
        %n_r = mod( (1:N_ag_r_b) + nodeShift - 1, msh.bandData.N_ag_r) + 1 + N_ag_s; %rotor nodes
        n_r = msh.bandData.sortedNodes_rotor( mod( (1:N_ag_r_b) + nodeShift - 1, msh.bandData.N_ag_r) + 1 );

        inds_n = [n_s n_r];
        ptemp = [msh.bandData.p_ag_virt(:, n_s) ...
            [cos(rotorAngle) -sin(rotorAngle); sin(rotorAngle) cos(rotorAngle)]*msh.bandData.p_ag_virt(:,n_r)];
        
        %triangulation
        ttemp = delaunay(ptemp(1,:), ptemp(2,:))';
        
        %removing elements inside rotor with inpolygon
        %pc = ( ptemp(:,ttemp(1,:)) + ptemp(:,ttemp(2,:)) + ptemp(:,ttemp(3,:)) ) / 3;
        %insidePolygon = msh.bandData.p_ag_virt(:, (N_ag_s+1):end);
        %Idx = inpolygon(pc(1,:), pc(2,:), insidePolygon(1,:), insidePolygon(2,:));
        %ttemp = ttemp(:, ~Idx);
        
        %removing elements consisting only of rotor nodes
        ttemp = ttemp(:, sum( ttemp <= numel(n_s) , 1)>0 );
        
        %removing "jumping" elements (that have nodes on both periodic
        %boundaries)
        el_1 = floor( (find(ttemp(:) == size(ptemp,2))-1) / 3) + 1; %elements touching 1st periodic boundary
        el_2 = floor( (find(ttemp(:) == (numel(n_s)+1))-1) / 3) + 1; %elements touching 2nd periodic boundary
        ttemp = ttemp(:, setdiff(1:size(ttemp,2), intersect(el_1, el_2)));
        
        %reverting to semi-global indexing (= corresponding to p_ag_virt)
        ttemp = reshape( inds_n(ttemp(:)), 3, []);

        %t_ag = msh.bandData.t_ag; t_ag(:, el_r) = ttemp;
        t_ag = [msh.bandData.t_ag(:, setdiff(1:msh.bandData.Ne_ag, el_r)) ttemp];
        p_ag = msh.bandData.p_ag_virt; p_ag(:, inds_n) = ptemp;
        
        S_ag_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', 1/mu0, [], struct('t', t_ag, 'p', p_ag), []);
        
        inds = 1:(S_ag_struct.ri-1);
        S_ag_struct.E(inds) = S_ag_struct.E(inds) .* transpose(msh.bandData.el_table(3, S_ag_struct.I(inds)));
        S_ag_struct.I(inds) = msh.bandData.el_table(2, S_ag_struct.I(inds));
        S_ag_struct.E(inds) = S_ag_struct.E(inds) .* transpose(msh.bandData.el_table(3, S_ag_struct.J(inds)));
        S_ag_struct.J(inds) = msh.bandData.el_table(2, S_ag_struct.J(inds));
        
        %assembling matrix
        S_ag = sparseFinalize(S_ag_struct, N, N);
        
    end
elseif numel(varargin) == 2
    error('Not usable with this implementation.');
end

end