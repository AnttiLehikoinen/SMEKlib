function sim = sim_setRotorCircuitMatrices(sim)
%sim_setRotorCircuitMatrices set rotor circuit.
%
% (c) 2017 Antti Lehikoinen / Aalto University

%setting FE-side circuit matrices
Mc = MatrixConstructor();
if sim.dims.type_rotorWinding == defs.cage
    %mass matrix for rotor
    rotorConductors = sim.msh.namedElements.get('rotorConductors');
    Mc.assemble_matrix(Nodal2D(Operators.I), Nodal2D(Operators.I), ...
        sim.dims.sigma_rotor, horzcat(rotorConductors{:}), sim.msh);
    %Mc = MatrixConstructor(Nodal2D(Operators.I), Nodal2D(Operators.I), ...
    %    sim.dims.sigma_rotor, horzcat(rotorConductors{:}), sim.msh);

    shaft = sim.msh.namedElements.get('shaft');
    if any(shaft)
        Mc.assemble_matrix(Nodal2D(Operators.I), Nodal2D(Operators.I), ...
           6e6, shaft, msh);
    end
    
    %rotor winding matrix
    Nc_r = numel(rotorConductors);
    %JF_struct = [];
    JF_s = MatrixConstructor;
    for k = 1:Nc_r
        %JF_struct = assemble_vector('', 'nodal', sim.dims.sigma_rotor, k, rotorConductors{k}, sim.msh, JF_struct);
        JF_s.assemble_vector(Nodal2D(Operators.I), k, sim.dims.sigma_rotor, rotorConductors{k}, sim.msh);
    end
    %sim.matrices.Cr = sparseFinalize(JF_struct, sim.Np, Nc_r);
    sim.matrices.Cr = JF_s.finalize(sim.Np, Nc_r);

    %conductor areas for rotor
    cAr = sum(sim.matrices.Cr*eye(Nc_r)/(sim.dims.sigma_rotor),1);
    sim.matrices.DRr = sparsediag(sim.dims.leff ./(sim.dims.sigma_rotor*cAr)); %resistance matrix
else
    %error('Invalid rotor winding type');
end

%checking if solid conductors elsewhere in the rotor
solids = sim.msh.namedElements.get('solidConductors_rotor');
if ~isempty( solids )
    for ks = 1:size(solids, 2)
        Mc.assemble_matrix(Nodal2D(Operators.I), Nodal2D(Operators.I), ...
           solids{1, ks}, solids{2, ks}, sim.msh);
    end
end

sim.matrices.Mr = Mc.finalize();

%setting loop matrices
if sim.dims.type_rotorWinding == defs.cage
    %loop matrix for rotor
    r_er = 1/sim.dims.Qr*pi*sim.dims.D_ring/(sim.dims.sigma_rotor*sim.dims.A_ring); %end-ring segment resistance
    [Lr, Zr] = rotorConnectionMatrix(sim.dims.Qr, sim.msh.symmetrySectors, sim.dims.p, 0, r_er);
    sim.matrices.Lr = Lr;
    sim.matrices.Zew_r = Zr;
else
    %error('Invalid rotor winding type');
end