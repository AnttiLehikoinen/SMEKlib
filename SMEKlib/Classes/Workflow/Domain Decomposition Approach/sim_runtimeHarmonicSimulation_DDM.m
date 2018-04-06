function sim = sim_runtimeHarmonicSimulation_DDM(sim, pars, varargin)
%sim_runtimeHarmonicSimulation_DDM domain-reduced time-harmonic simulation.
%
% (c) 2018 Antti Lehikoinen / Aalto University

f = pars.f;
%U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
U = 400/2
w = 2*pi*f;

if ~isempty(pars.slip)
    slips = pars.slip;
else
    slips = sim.dims.slip;
end
slip = slips(1);

%setting harmonic reluctivity function
if ~numel(varargin)
    %nu_struct = initialize_harmonicReluctivityStruct_interp1(sim.msh, true);
    nu_struct = initialize_reluctivityStruct_interp1(sim.msh, true);
    nu_fun = @(B)( calculate_reluctivity(B, nu_struct) );
else
    nu_fun = varargin{1};
end

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

%(rotor) circuit matrices
[Stot, Mtot] = get_circuitMatrices(sim, slip);

%reduced-basis coupling blocks
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;
P_m2D = sim.msh.misc.P_m2D;
L_s = sim.matrices.Ls;
tempcell_DD = cell(1, Qs_sector); [tempcell_DD{:}] = deal( sim.misc.R_AA );
tempcell_DI = cell(1, Qs_sector); [tempcell_DI{:}] = deal( sim.misc.R_AI );
tempcell_ID = cell(1, Qs_sector); [tempcell_ID{:}] = deal( sim.misc.R_UA );
tempcell_II = cell(1, Qs_sector); [tempcell_II{:}] = deal( sim.misc.R_UI );
Q_DD = P_m2D'*blkdiag(tempcell_DD{:})*P_m2D; clear tempcell_DD;
Q_DI = P_m2D'*blkdiag(tempcell_DI{:})*L_s; clear tempcell_DI;
Q_ID = L_s'*blkdiag(tempcell_ID{:})*P_m2D; clear tempcell_ID;
Q_II = L_s'*blkdiag(tempcell_II{:})*L_s; clear tempcell_II;

%numbers of variables
Nu = sim.results.Nu_r;
Ni_r = sim.results.Ni_r;
Ni_s = sim.results.Ni_s;

Stot = Stot + sim.msh.get_AGmatrix(0, size(Stot,1));
Qred = [Q_DD sparse(sim.Np, Nu) Q_DI sparse(sim.Np, Ni_r);
    sparse(Nu, sim.Np + Nu + Ni_r + Ni_s);
    Q_ID sparse(Ni_s,Nu) Q_II];

%temp solution
phi0 = 0;
UH = kron(eye(3), ones(size(L_s,2)/3,1)) * U*[exp(1i*phi0); exp(1i*phi0 - 1i*2*pi/3); exp(1i*phi0 - 1i*4*pi/3)];

Q = [Stot+real(Qred) w*Mtot+imag(Qred);
    -w*Mtot-imag(Qred) Stot+real(Qred)];
Nui = size(Stot,1) - size(sim.matrices.P,1);
PTT = blkdiag(sim.matrices.P, speye(Nui, Nui), ...
    sim.matrices.P, speye(Nui, Nui));

sim.misc.Q = Q;

%assembling load vector
Ftemp = [zeros(sim.Np + sim.results.Nu_r+sim.results.Nu_s,1); UH];
Ftot = [real(Ftemp); -imag(Ftemp)];

Xtot = zeros(size(Q,1), 1);
for kiter = 1:15
    % assembling Jacobian and residual blocks
    %[J11, J12, J21, J22, res11, res22] = assemble_ComplexJacobian(sim.nu_fun, Xtot(:,kslip), [], sim.msh);
    [J11, J12, J21, J22, res11, res22] = Jc.eval_complex(Xtot, nu_fun);
    
    sim.misc.res11 = res11;
    sim.misc.res22 = res22;
    sim.misc.J11 = J11;
    
    %finalizing
    Jtot = PTT'*( [J11 J12; J21 J22] + Q )*PTT;
    res_tot = PTT'*( Q*Xtot - Ftot + [res11;res22] );
    
    %checking convergence
    disp( norm(res_tot) / norm(Ftot) )
    if (norm(res_tot) / norm(Ftot)) < 1e-9
        break;
    end
    
    dX = - Jtot \ res_tot;
    Xtot = Xtot + PTT*dX;
end
sim.results.Xh = Xtot;

end