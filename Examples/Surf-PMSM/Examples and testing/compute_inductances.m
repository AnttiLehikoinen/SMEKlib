%angles = linspace(0, pi/2, 24);
IA = 10;
angles = linspace(0, pi, 21); phi0 = 0;
Imain = 10.5*[cos(angles + phi_bias + phi0); 
    cos(angles-2*pi/3 + phi_bias + phi0); 
    cos(angles - 4*pi/3 + phi_bias + phi0)];

%for testing
%angles = phi_bias*0; Imain = 12*[cos(angles); cos(angles-2*pi/3); cos(angles-4*pi/3)]*0;


%setting simulation parameters for a sweep
pars = SimulationParameters('Is', Imain, 'slip', 1, 'rotorAngle', zeros(1, size(Imain,2)));

dimsc2 = dimsc; dimsc2.supply_type = defs.current_supply;
sim_DC = MachineSimulation(mshc, dimsc2); %simulation object

sim_DC.run_static(pars);
A = sim_DC.results.Xs;

%computing forces and torque
T = sim_compute_torque(sim_DC, pars, 'static');

%T2 = compute_Torque_simple(A, mshc, 0.5*dim.Rout+0.5*dim.Sin, dim.Sin, [], pars.rotorAngle) * mshc.symmetrySectors*dimsc.leff;

figure(6); clf; hold on; box on;
plot(angles/pi*180, T, 'b');
mean(T)*2*pi*25

figure(5); clf; hold on; box on;
[~,ind] = min( abs(T - 10e3/(2*pi*25)) );
[~,inds] = sort( abs(T - 10e3/(2*pi*25)) );
%ind = 1
sim_DC.fluxplot(ind, pars, 'static');
%plot([0 dim.Sout], [0 dim.Sout], 'r');



%computing dq-axis stuff
ra = phi_bias;
rotM = [cos(-ra) -sin(-ra);
    sin(-ra) cos(-ra)];
angles_p = -[0 2/3 4/3]*pi;
Tp = 2/3*[cos(angles_p); -sin(angles_p)];
Phidq = rotM*Tp*sim_DC.matrices.Ls'*dimsc.leff*sim_DC.matrices.Cs'*A(1:sim_DC.Np,:) * mshc.symmetrySectors;
Idq = rotM*Tp*Imain;

figure(7); clf; hold on;
plot(Idq(1,:), Phidq(1,:), 'b' )
plot(Idq(2,:), Phidq(2,:), 'r' )
pd = polyfit(Idq(1,:), Phidq(1,:), 1);
pq = polyfit(Idq(2,:), Phidq(2,:), 1);
pd(1)/pq(1)
Ld = pd(1);
Lq = pq(1);