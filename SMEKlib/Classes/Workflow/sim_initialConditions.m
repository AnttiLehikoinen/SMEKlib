pars_ref = SimulationParameters('f', 1000, 'U', 400, 'slip', 0.5e-2, 'N_periods', 6);

pars = pars_ref;
sim = sim_ref;


f = pars.f;
U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
%U = 400/sqrt(3)/2 *sqrt(2) ; %taking the second-turn side outside symmetry sector into account
w = 2*pi*f;

if ~isempty(pars.slip)
    slips = pars.slip;
else
    slips = sim.dims.slip;
end

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

Ntot = size(sim.results.Xh, 1) / 2;
Np = sim.Np;
Nui = Ntot - size(sim.matrices.P,1);
PT = sim.matrices.P;

kslip = 1;

slip = slips(kslip);
[Stot, Mtot] = get_circuitMatrices(sim, slip);

Finit = -( Stot*sim.results.Xh(1:Ntot, kslip) + w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) );
Finit = Finit(1:Np,1);

X0 = sim.results.Xh(1:Np, kslip);
Sag = sim.msh.get_AGmatrix(0);
for kiter = 1:15
    [J, res] = Jc.eval(X0, sim.nu_fun);
    
    Jtot = PT'*( J + Sag )*PT;
    res_tot = PT'*(Sag*X0 - Finit + res );
    
    %checking convergence
    disp( norm(res_tot) / norm(Finit) )
    if (norm(res_tot) / norm(Finit)) < 1e-6
        break;
    end
    
    %Newton step
    dX = - Jtot \ res_tot;
    
    X0 = X0 + PT*dX;
end
X0 = [X0; sim.results.Xh((Np+1):Ntot, 1)];

return

%{
figure(5); clf;
%msh_trimesh(msh, sim.results.Xh(1:Np, kslip), []);
A = sim.results.Xh(1:Np, kslip);
drawFluxDensity(msh, A, 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 1.8]);
drawFluxLines(msh, A, 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);

figure(6); clf;
%msh_trimesh(msh, X, []);
A = X;
drawFluxDensity(msh, A, 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 1.8]);
drawFluxLines(msh, A, 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);
%}

sim.results.X0 = sim.results.Xh(1:Ntot,1); sim.results.X0(1:Np) = X0;

sim.run_timestepping(pars_ref);

figure(4); %clf; hold on; box on;
plot( pars_ref.ts, sim_ref.Is );
plot( pars_ref.ts, sum(sim_ref.Is,1), 'k--' );