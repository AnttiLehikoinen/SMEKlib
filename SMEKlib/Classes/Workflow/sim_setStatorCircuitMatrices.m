function sim = sim_setStatorCircuitMatrices(sim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stator winding matrices

if isfield(sim.dims, 'W')
    %winding configuration matrix given
    sim.matrices.W = sim.dims.W;
elseif isfield(sim.dims, 'N_layers')
    W = windingConfiguration_1(sim.dims.q, sim.dims.p, sim.dims.a, sim.dims.c);
    temp = unique( abs(W(:, 1:(sim.dims.Qs/sim.msh.symmetrySectors))) );
    if mod(numel(temp), sim.dims.a)
        warning('Symmetrizing winding for the symmetry sector');
        W = (floor( (abs(W)-1)/sim.dims.a )+1).*sign(W);
    end
    sim.matrices.W = W;
else
    sim.matrices.W = windingConfiguration_1(sim.dims.q, sim.dims.p);
end

%JF_struct = [];
JF_c = MatrixConstructor;
statorConductors = sim.msh.namedElements.get('statorConductors');
if sim.dims.type_statorWinding == defs.stranded
    Nc_s = numel(statorConductors);
    for k = 1:Nc_s
        %JF_struct = assemble_vector('', 'nodal', 1, k, statorConductors{k}, sim.msh, JF_struct);
        JF_c.assemble_vector(Nodal2D(Operators.I), k, 1, statorConductors{k}, sim.msh);
    end
    %JF_s = sparseFinalize(JF_struct, sim.Np, Nc_s);
    JF_s = JF_c.finalize(sim.Np, Nc_s);
    cAs = sum(JF_s*eye(Nc_s),1);
    if isfield(sim.dims, 'l_halfCoil')
        DRs = sparsediag(sim.dims.l_halfCoil ./(sim.dims.sigma_stator*cAs)) / sim.dims.fillingFactor; %resistance matrix
    else
        DRs = sparsediag(sim.dims.leff ./(sim.dims.sigma_stator*cAs)) / sim.dims.fillingFactor; %resistance matrix
    end
    sim.matrices.DRs = DRs;
    sim.matrices.Cs = bsxfun(@times, JF_s, 1./cAs);    
    sim.matrices.Ms = sparse(sim.Np, sim.Np);    
    
    %loop matrix
    Ls = statorConnectionMatrix(sim.matrices.W, size(sim.matrices.W,1), 1);
    Ls = Ls(1:(size(Ls,1)/sim.msh.symmetrySectors),:); 
    sim.matrices.Ls = Ls(:, sum(abs(Ls),1)>0) * sim.dims.N_series;
    sim.matrices.Zew_s = sparse(sim.matrices.Ls'*DRs*sim.matrices.Ls);
elseif sim.dims.type_statorWinding == defs.decomposed
    % decomposed winding: A-parts handled elsewhere
     %loop matrix
    Ls = statorConnectionMatrix(sim.matrices.W, sim.dims.Nc_slot, sim.dims.N_series, 4, 3);
    %Ls = statorConnectionMatrix(sim.matrices.W, sim.dims.Nc_slot, sim.dims.N_series);
    Ls = Ls(1:(size(Ls,1)/sim.msh.symmetrySectors),:); 
    sim.matrices.Ls = sparse(Ls(:, sum(abs(Ls),1)>0));
    sim.matrices.Ms = sparse(sim.Np,sim.Np);
    sim.matrices.Cs = sparse(sim.Np, size(Ls,1));
else
    warning('Stator winding type not yet implemented.');
end

%end-winding inductance given?
if isfield(sim.dims, 'Lew')
    Lew = sim.dims.Lew;
    if numel(Lew) == 1
        sim.matrices.Zew_s = sim.matrices.Zew_s + 1i*Lew*eye(size(sim.matrices.Zew_s,1));
    else
        sim.matrices.Zew_s = sim.matrices.Zew_s + 1i*Lew;
    end
end

end