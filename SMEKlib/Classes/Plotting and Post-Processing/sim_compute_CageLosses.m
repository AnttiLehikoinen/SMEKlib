function [Pmean, P_bar, J] = sim_compute_CageLosses(sim, pars, varargin)
%fun_PMlosses Compute PM losses from results.
%
% Call syntax
% [Pmean, P, J] = sim_compute_CageLosses(sim, pars)
% [Pmean, P, J] = sim_compute_CageLosses(sim, pars, plotting)
% [Pmean, P, J] = sim_compute_CageLosses(sim, pars, plotting, steps)
%
% where
%   Pmean = losses (W) in each bar, at each time-step.
%   J = current density (A/m^2) at each node at each time-step. Correct
%       only when different PMs don't share nodes.
%
% and
%   sim = MachineSimulation object
%   pars = SimulationParameters object
%   plotting = boolean, plot losses at each time-step? (Optional)
%   steps = integer array, plot losses at these steps only? Only active
%       when plotting = true. (Optional)
%
% Note: breaks down if there are other conductors in the rotor, for now.
%
% (c) 2018-2019 Antti Lehikoinen / Smeklab

MAKE_ANIMATION = false; %save J-animation into a gif file?

%plot losses?
if numel(varargin)
    plotting_on = varargin{1};
else
    plotting_on = false;
end

ts = pars.ts;
Nsamples = numel(ts);
Np = sim.Np; Nu_s = sim.results.Nu_s;
msh = sim.msh;

bars = sim.msh.namedElements.get('rotorConductors');

n_bars = cell(1, numel(bars));
M_bars = cell( size(n_bars) );
for k = 1:numel(n_bars)
    n_bars{k} = toRow( unique( msh.t(:, bars{k}) ) );
    M_bars{k} = MatrixConstructor(Nodal2D, Nodal2D, sim.dims.leff, bars{k}, msh).finalize();
end
sigma_rotor = sim.dims.sigma_rotor;
if numel(sigma_rotor) == 1
    sigma_rotor = sigma_rotor*ones(1, numel(n_bars));
end

alpha2 = pars.alpha2; alpha1 = 2.0 - pars.alpha2;
P_bar = zeros(2, Nsamples);
J = zeros(Np, Nsamples);


%for animation
rotorAngles = 2*pi*pars.f*ts / sim.dims.p;
for ks = 2:Nsamples
    %disp(num2str(ks));
    
    dt = ts(ks) - ts(ks-1);
    
    Apm = sim.results.Xt(:,ks);
    dApm = (Apm - sim.results.Xt(:,ks-1) ) / dt;
    
    
    for kpm = 1:numel(n_bars)
        Jpm = zeros(Np, 1);
        Jpm(n_bars{kpm}) = -dApm( n_bars{kpm} ) * sigma_rotor(kpm);
        
        Ueff = ( alpha2*sim.results.Xt(Np+Nu_s+kpm, ks) + ...
            alpha1*sim.results.Xt(Np+Nu_s+kpm, ks-1) ) / 2;
        
        Jpm(n_bars{kpm}) = Jpm(n_bars{kpm}) +  sigma_rotor(kpm)*Ueff/sim.dims.leff;
        
        P_bar(kpm, ks) = Jpm' * M_bars{kpm} * (Jpm/sigma_rotor(kpm));
        J(:, ks) = J(:, ks) + Jpm;
    end
end
P_bar = P_bar * sim.msh.symmetrySectors;

Pmean = mean(P_bar, 2);

if plotting_on
    if numel(varargin) == 2
        ks_vec = varargin{2};
    else
        ks_vec = 2:Nsamples;
    end
    
    filename = 'cage_currents.gif'; ri = 1;
    h = figure(10); clf; hold on; box on;
    colorbar;
    cax = 1e-6*[min(min(J(:, ks_vec))) max(max(J(:, ks_vec)))];
    
    %getting axis limits
    if sim.msh.symmetrySectors == 4
        axlims = [-sim.dims.Rout sim.dims.Sout -sim.dims.Rout sim.dims.Sout];
    else
        axlims = sim.dims.Sout*[-1 1 -1 1];
    end
    
    Xplot = zeros(size(msh.t,1), size(msh.t,2));
    Yplot = Xplot;
    for kp = 1:size(Xplot, 1)
        Xplot(kp,:) = msh.p(1, msh.t(kp,:));
        Yplot(kp,:) = msh.p(2, msh.t(kp,:));
    end
    if size(Xplot,1)==6
        ord = [1 4 2 5 3 6];
    else
        ord = 1:3;
    end
    
    for ks = ks_vec
        clf;
        msh_triplot(msh, [], rotorAngles(ks), 'k');
        drawCurrentDensity(msh, 1e-6*J(:,ks), [bars{:}], rotorAngles(ks), 'linestyle', 'none');
        colorbar; caxis(cax);
        title('Current density (Apeak / mm^2)');
        axis(axlims);
        drawnow;
        
        if MAKE_ANIMATION
            %axis([-0.07 0.1 -0.07 0.08]); daspect([1 1 1]);
            %axis equal;
            %drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            if ri == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 1e-3);
            elseif ri < numel(ks_vec)
                imwrite(imind,cm,filename,'gif', 'WriteMode','append', 'DelayTime', 1e-3);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1.5);
            end
            ri = ri + 1;
        end
    end
end

end