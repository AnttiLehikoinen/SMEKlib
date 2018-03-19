function [Bx, By, w_quad, DETF] = sim_getQuadB(sim)
%sim_getB returns flux density at integration points.
%
% [Bx, By, w_quad, DETF] = sim_getQuadB(sim)
% returns Bx and By at the integration points, with the corresponding
% quadrature weights w_quad and mapping determinants DETF, corresponding to
% the time-stepping data contained in MachineSimulation object sim.
%
% Only works for time-stepping results.
%
% Bx and By are 3D-arrays, with the following indexing
% Bx(k_quad, k_elements, k_timestep)
%
% (c) 2018 Antti Lehikoinen / Aalto University


%Jacobian constructor for quick access to shape function values
Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), true);

N_shape = size(Jc.Fvals_shape, 2);
Ne = size(sim.msh.t,2);

DETF = Jc.DETF;
	
Nsamples = size(sim.results.Xt, 2);
N_quad = numel(Jc.w);
w_quad = Jc.w;

Bx = zeros(N_quad, Ne, Nsamples);
By = zeros(N_quad, Ne, Nsamples);
X = sim.results.Xt;
for k_shape = 1:N_shape
	X_shape = reshape( X(sim.msh.t(k_shape,:), :), 1, Ne, [] );
	for k_quad = 1:N_quad
		Bx(k_quad, :, :) = Bx(k_quad, :, :) + bsxfun(@times, Jc.Fvals_shape{k_quad, k_shape}(1,:), X_shape );
		By(k_quad, :, :) = By(k_quad, :, :) + bsxfun(@times, Jc.Fvals_shape{k_quad, k_shape}(2,:), X_shape );
	end
end

end