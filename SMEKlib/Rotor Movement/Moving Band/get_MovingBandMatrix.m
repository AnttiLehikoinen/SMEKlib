function S_ag = get_MovingBandMatrix(rotorAngle, varargin)
%get_MovingBandMatrix moving band stiffness matrix.
% 
% S_ag = get_MovingBandMatrix(rotorAngle, msh) returns the moving band
% matrix of size Np x Np.
% Note: all moving band elements are shifted at the same time.
% 
% S_ag = get_MovingBandMatrix(rotorAngle, msh, N) returns a zero-appended
% larger matrix of size NxN, with N>=Np.
%
% Copyright (c) 2015-2016 Antti Lehikoinen / Aalto University


mu0 = pi*4e-7;
if (numel(varargin) == 1) || ( (numel(varargin)==2) && isa(varargin{2}, 'double'))
    msh = varargin{1};
    
    if ( (numel(varargin)==2) && isa(varargin{2}, 'double'))
        N = varargin{2};
    else
        N = size(msh.p,2);
    end
    
    if ~isfield(msh, 'symmetrySectors')
        %full air-gap mesh
        
        %special rotation function supplied?
        if isfield(msh, 'updateRotorPosition')
            [t_ag, ~, pnew] = msh.updateRotorPosition(rotorAngle);
        else
            [t_ag, ~, pnew] = updateRotorPosition(rotorAngle, msh);
        end

        msh_temp = struct('t', t_ag, 'p', pnew);

        S_ag_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', 1/mu0, [], msh_temp, []);
        S_ag = sparseFinalize(S_ag_struct, size(msh.p,2), size(msh.p,2));
    else
        %number of nodes to skip on rotor surface
        nodeShift = -floor( (rotorAngle - 1*msh.bandData.shiftTol/2) / msh.bandData.shiftTol ) - 1;

        %shifting indices in the sorted list
        newPositions = mod(msh.bandData.originalPositions_rotor - 1 + nodeShift, msh.bandData.N_ag_r ) + 1;

        t_ag = msh.bandData.t_ag;
        t_ag(msh.bandData.inds_r) = msh.bandData.sortedNodes_rotor(newPositions);
        
        p = msh.bandData.p_ag_virt;
        p(:, msh.bandData.sortedNodes_rotor) = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)] * ...
            p(:, msh.bandData.sortedNodes_rotor);
        
        %figure(3); clf;
        %triplot(t_ag', p(1,:), p(2,:));
        
        S_ag_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', 1/mu0, [], struct('t', t_ag, 'p', p), []);
        
        %figure(11); clf; msh_triplot(struct('t', t_ag, 'p', p), [], 'b');
        
        % handling periodicity conditions
        %S_ag_struct = sparseEliminate(S_ag_struct, true, true, msh.bandData.el_table)        
        inds = 1:(S_ag_struct.ri-1);
        S_ag_struct.E(inds) = S_ag_struct.E(inds) .* transpose(msh.bandData.el_table(3, S_ag_struct.I(inds)));
        S_ag_struct.I(inds) = msh.bandData.el_table(2, S_ag_struct.I(inds));
        S_ag_struct.E(inds) = S_ag_struct.E(inds) .* transpose(msh.bandData.el_table(3, S_ag_struct.J(inds)));
        S_ag_struct.J(inds) = msh.bandData.el_table(2, S_ag_struct.J(inds));
        
        %assembling matrix
        S_ag = sparseFinalize(S_ag_struct, N, N);        
    end
elseif numel(varargin) == 2
    % polynomial approximation of the matrix; not yet in use.
    
    msh = varargin{1};
    bandElementContributions = varargin{2};
    
    Ne_ag = msh.bandData.Ne_ag;
    
    t_ag = msh.updateRotorPosition(rotorAngle);
    
    %generating references (1st order mesh & functions assumed)
    temp = repmat(1:3,3,1);
    row_inds_ref = reshape(temp',1,[]);
    column_inds_ref = temp(:)';
    
    %calculating air-gap node angles
    inds_rotor = (msh.bandData.N_ag_s + 1):(msh.bandData.N_ag_s + msh.bandData.N_ag_r);
    agAngles = msh.bandData.agAngles_all;
    agAngles(:, inds_rotor) = agAngles(:, inds_rotor) + rotorAngle;
    
    %calculating pu element angles
    angle_element = angleDifference( agAngles(:, t_ag(1,:) ), agAngles(:, t_ag(3,:) ) )' ;
    angle_element_pu = bsxfun(@plus, angle_element, - bandElementContributions.element_a0);
    angle_element_pu = bsxfun(@times, angle_element_pu, 1./bandElementContributions.element_ad);
    
    %element determinant values
    pdet = aux_polyvals(bandElementContributions.elementPolynomials_det, angle_element_pu);
    
    I_ag = zeros(9*Ne_ag, 1); J_ag = I_ag; E_ag = I_ag;
    
    for k_mp = 1:9
        %computing element contributions
        elementContributions = aux_polyvals(bandElementContributions.elementPolynomials_entr{k_mp}, angle_element_pu) ./ pdet;
        
        %storing them for as sparse-triplets
        inds = ((k_mp-1)*Ne_ag + 1):(k_mp*Ne_ag);        
        E_ag(inds) = elementContributions;
        I_ag(inds) = t_ag(row_inds_ref(k_mp), :)';
        J_ag(inds) = t_ag(column_inds_ref(k_mp), :)';
    end
    
    S_ag = sparse(msh.bandData.agNodes_global(I_ag), msh.bandData.agNodes_global(J_ag), E_ag/mu0, size(msh.p,2), size(msh.p,2));    
end


end