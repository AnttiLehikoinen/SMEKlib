%d,q current plots

ts = pars.ts;

rotorAngles = wm*ts;
Is = sim.Is;

Idq = zeros(2, numel(ts));
for k = 1:numel(ts)
    ra = phi_bias + rotorAngles(k)*dimsc.p;
    rotM = [cos(-ra) -sin(-ra);
        sin(-ra) cos(-ra)];
    angles_p = -[0 2/3 4/3]*pi;
    Tp = 2/3*[cos(angles_p); -sin(angles_p)];
    %Phidq = rotM*Tp*sim.matrices.Ls'*sim.matrices.Cs'*A(1:simc.Np,:);
    Idq(:,k) = rotM*Tp*Is(:,k);
end

figure(12); clf; hold on;
plot(ts, Idq(1,:), 'b');
plot(ts, Idq(2,:), 'r');

figure(13); clf; hold on;
Uhere = 520;
phi0 = phi_bias-pi/2 - pi/180*15;
Uplot = Uhere/sim.msh.symmetrySectors * sim.dims.a * sqrt(2)* ...
    [cos(w*ts-phi0); cos(w*ts - 2*pi/3-phi0); cos(w*ts - 4*pi/3-phi0)];

plot(Uplot(1,:), 'b');
plot(Uplot(2,:), 'r');
plot(Uplot(3,:), 'k');

plot(ctrl.Us(1,:), 'b');
plot(ctrl.Us(2,:), 'r');
plot(ctrl.Us(3,:), 'k');