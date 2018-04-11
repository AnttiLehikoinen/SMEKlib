Xh = sim.results.Xh; %harmonic solution

%back to complex domain
Xh = reshape(Xh, [], 2); Xh = Xh(:,1) - 1i*Xh(:,2);

Is = Xh( (end-2):end, :);
Ur = Xh(sim.Np+(1:sim.results.Nu_r));
Ir = sim.matrices.Zew_r \ (-sim.matrices.Lr'*Ur);

figure(1); clf;
compass(real(Is), imag(Is))

figure(2); clf;
compass(real(Ir), imag(Ir))

%linearized reluctivity
[~, B] = calculate_B(real(Xh), mshc);
nu_real = calculate_reluctivity(B, sim.nu_struct);

[~, B] = calculate_B(imag(Xh), mshc);
nu_imag = calculate_reluctivity(B, sim.nu_struct);
nu_lin = sqrt(nu_real.^2 + nu_imag.^2);

S = MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), nu_lin, [], mshc).finalize(Np, Np);

%voltages
w = 2*pi*f;
slip = dims.slip;
M = 0*1i*w*(sim.matrices.Ms + slip*sim.matrices.Mr);

Cr = sim.matrices.Cr;
cAr = sum(Cr*eye(size(Cr,2)),1);
Cr = bsxfun(@times, Cr, 1./cAr);


A = (sim.matrices.P'*(S + M + mshc.get_AGmatrix(0))*sim.matrices.P) \ ...
    (sim.matrices.P'*[sim.matrices.Cs*sim.matrices.Ls ...
    Cr*sim.matrices.Lr]);
A = sim.matrices.P*A;

%{
figure(3); clf; hold on;
k = 4;
drawFluxDensity(mshc, A(:,k), 'LineStyle', 'none'); colormap('jet'); colorbar; %caxis([0 2]);
drawFluxLines(mshc, A(:,k), 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);
%}


U = [1i*w*sim.dims.leff*sim.matrices.Ls'*sim.matrices.Cs';
    1i*slip*w*sim.dims.leff*sim.matrices.Lr'*Cr']*A + blkdiag(0*sim.matrices.Zew_s, 0*sim.matrices.DRr*sim.matrices.Lr);
%U = imag(U);

%return

Ps = 1/3*[1 exp(1i*2*pi/3) exp(1i*4*pi/3)];
PsI = [1 exp(1i*2*pi/3) exp(1i*4*pi/3)]';
Pr = exp(1i*(0:9)*pi/10) * 1/10 *4;
PrI = exp(1i*(0:9)*pi/10)' *4;

%PsI = ones(3,1); PrI = ones(10,1);

Ppark = blkdiag(Ps, Pr);
PparkI = blkdiag(PsI, PrI);

%figure(3); clf;
%compass(real(Pr), imag(Pr))

%c1 = Pr*U(4:end,1:3)*PsI
%c2 = Ps*U(1:3,4:end)*PrI

Zp = Ppark*U*PparkI
Zp = Zp*diag([1 exp(1i*(pi/2-angle(Zp(1,2))))])
Zp = diag([1 Zp(1,2)/Zp(2,1)])*Zp
%scale = sqrt(Zp(2,2)/Zp(1,1));

%Zp = diag([1 scale])*Zp*diag([1 scale])