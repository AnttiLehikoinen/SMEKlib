%parse_mesh Combine stator and rotor mesh data into a mesh object.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parsing dimsc

dim.D_so = 2*dim.Sout; %needed for plotting
dim.q = dim.Qs/(3*2*dim.p); %number of slots per pole and phase

if ~isfield(dim, 'c')
    dim.c = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ne_r = size(datr.t_r, 2);
Np_r = size(datr.p_r, 2);
matel = [datr.m_r dats.m_s];
mshc = MachineMesh([datr.p_r dats.p_s], [datr.t_r Np_r+dats.t_s], matel);
mshc.setSymmetrySectors(dim.symmetrySectors, dim);

%setting special nodes
mshc.namedNodes.set( 'Dirichlet', [datr.n_dir_r Np_r+dats.n_dir_s] );
mshc.namedNodes.set( 'n_ag_s', Np_r+dats.n_ag_s );
mshc.namedNodes.set( 'n_ag_r', datr.n_ag_r );
mshc.namedNodes.set('Periodic_master', [datr.n_cw_r Np_r+dats.n_cw_s]);
mshc.namedNodes.set('Periodic_slave', [datr.n_ccw_r Np_r+dats.n_ccw_s]);

%setting special elements,
mshc.namedElements.set('rotel', 1:Ne_r); %rotor elements
rotorConductors = datr.RC;
mshc.namedElements.set('rotorConductors', rotorConductors);
statorConductors = cellfun(@(x)(x+Ne_r), dats.SC, 'UniformOutput', false);
mshc.namedElements.set('statorConductors', statorConductors);
%mshc.namedElements.set('PMs', datr.PMs);

%generating airgap triangulation and movement stuff
mshc.generateMovingBand(); %airgap triangulation

figure(1); clf; hold on; box on; axis equal tight;
ls = '-';
msh_triplot(mshc, [], 'b');
msh_fill(mshc, find(mshc.matel == dim.SM), [1 1 1]*0.7, 'linestyle', ls);
msh_fill(mshc, find(mshc.matel == dim.RM), [1 1 1]*0.7, 'linestyle', ls);

%msh_plot(mshc, mshc.namedNodes.get('Periodic_master'), 'gv-');
%msh_plot(mshc, mshc.namedNodes.get('Periodic_slave'), 'ro-');
%msh_plot(mshc, mshc.namedNodes.get('Dirichlet'), 'ko');
%msh_plot(mshc, mshc.namedNodes.get('n_ag_s'), 'cd-');
%msh_plot(mshc, mshc.namedNodes.get('n_ag_r'), 'md-');

%msh_fill(mshc, [statorConductors{:}], 'b', 'linestyle', 'none');
msh_fill(mshc, [rotorConductors{:}], 'y', 'linestyle', ls);

%plotting airgap triangulation
AGT = mshc.bandData;
[~, pag, tag] = mshc.bandData.t_ag(0);
triplot(tag', mshc.p(1,:), mshc.p(2,:), 'm');
%msh_triplot(AGT.msh_ag, [], 'm'); 

%default winding type plot
if dim.N_layers > 1
    Wtemp = windingConfiguration_1(dim.q, dim.p, 1, dim.c);
else
    Wtemp = windingConfiguration_1(dim.q, dim.p);
end
%phaseColors = {[139 69 19]/256, [0 0 0], [1 1 1]*0.6};
phaseColors = {[0 0 1], [1 0 0], [0 1 0]};
for kphase = 1:3
    inds_plus = find(Wtemp == kphase);
    inds_plus = inds_plus(inds_plus <= numel(statorConductors));
    if any(inds_plus)
        msh_fill(mshc, [statorConductors{inds_plus}], phaseColors{kphase}, 'linestyle', ls);
    end
    
    inds_minus = find(Wtemp == -kphase);
    inds_minus = inds_minus(inds_minus <= numel(statorConductors));
    if any(inds_minus)
        els = [statorConductors{inds_minus}];
        msh_fill(mshc, els(1:2:end), phaseColors{kphase}, 'linestyle', ls);
        msh_fill(mshc, els(2:2:end), phaseColors{kphase}*0.8, 'linestyle', ls);
    end
end

if true
    msh_plot(mshc, mshc.namedNodes.get('Periodic_master'), 'gv-');
    msh_plot(mshc, mshc.namedNodes.get('Periodic_slave'), 'ro-');
    msh_plot(mshc, mshc.namedNodes.get('Dirichlet'), 'ko');
    msh_plot(mshc, mshc.namedNodes.get('n_ag_s'), 'cd-');
    msh_plot(mshc, mshc.namedNodes.get('n_ag_r'), 'md-');
end