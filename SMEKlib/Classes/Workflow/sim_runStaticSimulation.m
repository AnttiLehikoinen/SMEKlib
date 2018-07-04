function sim = sim_runStaticSimulation(sim, pars)
%sim_runStaticSimulation general static simulation.
%
% WIP
%
% (c) 2018 Antti Lehikoinen / Smeklab


[Sc, ~] = get_circuitMatrices(sim);

Ntot = size(Sc, 1);
Nui = Ntot - size(sim.matrices.P,1);
PT = blkdiag(sim.matrices.P, speye(Nui));

%Jacobian constructor object
Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

S_ag = sim.msh.get_AGmatrix(0, Ntot);
Qconst = S_ag + Sc;

N = size(pars.U, 2);

FL = bsxfun(@plus, pars.U, sim.matrices.F);
Xs_all = zeros(Ntot, N);

for k = 1:N
    if k == 1
        Xs = zeros(Ntot, 1);
    else
        Xs = Xs_all(:,k-1);
    end
    for kiter = 1:50
        [J, res] = Jc.eval(Xs, sim.nu_fun);

        res_tot = PT'*(res + Qconst*Xs - FL(:,k));

        resNorm = norm(res_tot) / norm(FL(:,k));
        disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        %disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']); fflush(stdout); uncomment in Octave
        if resNorm < 1e-6
            break;
        end

        Jtot = PT' * (J + Qconst) * PT;

        Xs = Xs - PT*( Jtot\res_tot );
    end
    Xs_all(:,k) = Xs;
end

sim.results.Xs = Xs_all;

end