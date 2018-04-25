function sim = sim_runTimeHarmonicSimulation(sim, pars, varargin)
%run_timeHarmonicSimulation
%
% (c) 2017 Antti Lehikoinen / Aalto University


f = pars.f;
U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
%U = 400/sqrt(3)/2 *sqrt(2) ; %taking the second-turn side outside symmetry sector into account
w = 2*pi*f;

if ~isempty(pars.slip)
    slips = pars.slip;
else
    slips = sim.dims.slip;
end

%setting harmonic reluctivity function
if ~numel(varargin)
    nu_struct = initialize_harmonicReluctivityStruct_interp1(sim.msh, true);
    %nu_struct = initialize_reluctivityStruct_interp1(sim.msh, true);
    nu_fun = @(B)( calculate_reluctivity(B, nu_struct) );
else
    nu_fun = varargin{1};
end

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);
for kslip = 1:numel(slips)
    slip = slips(kslip);
    
    [Stot, Mtot] = get_circuitMatrices(sim, slip);
    Stot = Stot + sim.msh.get_AGmatrix(0, size(Stot,1));

    Q = [Stot w*Mtot;
        -w*Mtot Stot];

    %sim.matrices.Stot = Stot;
    %sim.matrices.Mtot = Mtot;
    %size(Stot)

    Nui = size(Stot,1) - size(sim.matrices.P,1);
    Nu = sim.results.Nu_r+sim.results.Nu_s;
    PTT = blkdiag(sim.matrices.P, speye(Nui, Nui), ...
        sim.matrices.P, speye(Nui, Nui));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %voltage vector
    t0 = 0;
    if numel(pars.U) > 1
        FI = pars.U;
    elseif sim.dims.connection_stator == defs.delta
        FI = U.*[exp(1i*w*t0); exp(1i*w*t0-1i*2*pi/3); exp(1i*w*t0-1i*4*pi/3)];
    else
        FI = U.* [exp(1i*w*t0); exp(1i*w*t0-1i*pi/3)];
    end

    %assembling load vector
    %Ftemp = [sim.matrices.F; zeros(Nu,1); FI(1:sim.results.Ni_s)];
    %Ftot = [real(Ftemp); -imag(Ftemp)];
    Ftot = [sim.matrices.F; zeros(Nu,1); real(FI(1:sim.results.Ni_s));
        -0*sim.matrices.F; zeros(Nu,1); -imag(FI(1:sim.results.Ni_s))];

    if kslip == 1
        Xtot = zeros(size(Q,1), numel(slips));
    end
    for kiter = 1:15
        % assembling Jacobian and residual blocks
        %[J11, J12, J21, J22, res11, res22] = assemble_ComplexJacobian(sim.nu_fun, Xtot(:,kslip), [], sim.msh);    
        [J11, J12, J21, J22, res11, res22] = Jc.eval_complex(Xtot(:,kslip), nu_fun); 

        %finalizing
        Jtot = PTT'*( [J11 J12; J21 J22] + Q )*PTT;    
        res_tot = PTT'*( Q*Xtot(:,kslip) - Ftot + [res11;res22] );

        %checking convergence
        disp( norm(res_tot) / norm(Ftot) )
        if (norm(res_tot) / norm(Ftot)) < 1e-6
            break;
        end

        dX = - Jtot \ res_tot;
        Xtot(:,kslip) = Xtot(:,kslip) + PTT*dX;
    end
    if kslip < numel(slips)
        Xtot(:,kslip+1) = Xtot(:,kslip);
    end
end
sim.results.Xh = Xtot;

end