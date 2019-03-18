%computing dq-axis stuff
ra = phi_bias;
rotM = [cos(-ra) -sin(-ra);
    sin(-ra) cos(-ra)];
angles_p = -[0 2/3 4/3]*pi;
Tp = 2/3*[cos(angles_p); -sin(angles_p)];
%Phidq = rotM*Tp*simc.matrices.Ls'*dimsc.leff*simc.matrices.Cs'*A(1:simc.Np,:) * mshc.symmetrySectors;
%Idq = rotM*Tp*Imain;

%rotorAngles = wm*pars.ts;
rotorAngles = pars.rotorAngle;
Nsamples = numel(rotorAngles);

Phidq = zeros(2, Nsamples);
Idq = zeros(2, Nsamples);

A = sim.results.Xs;
%Imain = sim.Is;
for k = 1:Nsamples
    ra = phi_bias + rotorAngles(k)*dimsc.p;
    rotM = [cos(-ra) -sin(-ra);
        sin(-ra) cos(-ra)];
    
    Phidq(:,k) = rotM*Tp*sim.matrices.Ls'*dimsc.leff*sim.matrices.Cs'*A(1:sim.Np,k) * mshc.symmetrySectors;
    Idq(:,k) = rotM*Tp*Imain(:,k);
end

figure(7); clf; hold on;
plot(rotorAngles, Phidq' );
%plot(rotorAngles, Phi_d(:, 56), 'bx-');
%plot(rotorAngles, Phi_q(:, 56), 'rx-');

figure(8); clf;
plot( Idq')

%computing
Id = Idq(1,:);
Iq = Idq(2,:);
Phi_d = Phidq(1,:);
Phi_q = Phidq(2,:);
c = [[ones(numel(Id),1) Id' 0*Iq']; 
    [zeros(numel(Id),1) 0*Id' Iq']] \ [Phi_d'; Phi_q']
%Tb = 3/2*dimsc.p * ( -Idq(1,:).*Phidq(2,:) + Idq(2,:).*Phidq(1,:) );
%figure(6); clf; hold on; box on; grid on;
%plot(angles/pi*180, T, 'b')
%plot(angles/pi*180, Tb, 'r')
    