

AGT = mshcs.bandData;

kt = 150;
Xref = sims.results.Xt(:, kt);
rA = rotorAngles(kt);

nmid = AGT.n_int;
angles = atan2(mshc.p(2,nmid), mshc.p(1,nmid));

%Xref(nmid) = sin(angles*2) + 0.05*cos(angles*50);

figure(7); clf; hold on; box on; grid on;
plot(angles/pi*180, Xref(nmid), 'bx-')



[Sint, P] = get_interpolatedMatrix(AGT, rA);

X2 = P*Xref(1:size(mshc.p,2));

%plot( mod(angles+rA, pi/4)/pi*180, X2(nmid), 'ro');
plot( mod(angles+rA, pi/4)/pi*180, X2(nmid)'.*( (-1).^(floor((angles+rA)/(pi/4))) ), 'ro-');

legend('Static', 'Interpolated');

xlabel('Angular coordinate (deg)');
ylabel('Vector potential (Wb/m)');