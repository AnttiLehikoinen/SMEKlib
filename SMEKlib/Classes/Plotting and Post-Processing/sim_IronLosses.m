function [Ptot, Physt, Peddy, Pexcess] = sim_IronLosses(sim, pars, varargin)
%fun_IronLosses Compute iron losses.
%
% Call syntax:
% [Ptot, Physt, Peddy, Pexcess] = fun_IronLosses(sim, pars)
% [Ptot, Physt, Peddy, Pexcess] = fun_IronLosses(sim, pars, plotting)
% [Ptot, Physt, Peddy, Pexcess] = fun_IronLosses(sim, pars, plotting, pretty_plotting)
%
% where
%   Ptot = total average iron losses
%   Physt, Peddy, Pexcess = element-wise average losses (W/kg)
%
% and
%   sim = MachineSimulation object
%   pars = SimulationParameters object
%   plotting = boolean, plot losses? (Optional)
%   pretty_plotting = boolean, generate smoother loss plot. Only active
%       with second-order elements and plotting=true. (Optional)
%
% (c) 2018 Antti Lehikoinen / Smeklab Ltd


if numel(varargin)
    plotting_on = varargin{1};
    if numel(varargin) == 2
        pretty_plotting = varargin{2};
    else
        pretty_plotting = false;
    end
else
    plotting_on = false;
    pretty_plotting = false;
end

if isempty(sim.results.Xt)
    A = sim.results.Xs(1:sim.Np,:);
    tsamples = pars.rotorAngle / (2*pi*pars.f / sim.dims.p);
else
    A = sim.results.Xt(1:sim.Np, :);
    tsamples = pars.ts;
end
%A(:,1) = A(:,2);

%not including all periods
%{
inds_t = 100:size(A, 2);
A = A(:, inds_t);
tsamples = tsamples(inds_t);
%}

msh = sim.msh;

fref = 50;

%preparing for fft
L = numel(tsamples);
Nsamples = L;
Fs = 1/( tsamples(2)-tsamples(1) );
Nrel = floor( 0.8*numel(tsamples)/2 );
%Nrel = 19;
freqs = Fs*(0:Nrel)/numel(tsamples);
freqRatios = freqs / fref;

Ncurl = Nodal2D(Operators.curl, Operators.curl);
Jc = JacobianConstructor(sim.msh, Ncurl, Ncurl, true);
N_quad = numel(Jc.w);
N_shape = size(Jc.Fvals_shape, 2);
Ne = size(sim.msh.t,2);

CB2 = cell(N_quad, 1);
for k_quad = 1:N_quad
    Bx = zeros(Ne, Nsamples);
    By = zeros(Ne, Nsamples);
    for k_shape = 1:N_shape
        X_shape = A(sim.msh.t(k_shape,:), :);
        
        Bx = Bx + bsxfun(@times, Jc.Fvals_shape{k_quad, k_shape}(1,:)', X_shape );
        By = By + bsxfun(@times, Jc.Fvals_shape{k_quad, k_shape}(2,:)', X_shape );
    end
    %fft
    Cx = fft(Bx, [], 2)/L*2;
    Cy = fft(By, [], 2)/L*2;
    CB2{k_quad} = abs(Cx).^2 + abs(Cy).^2;
end

%loss densities at nodal points?
if pretty_plotting && (msh.elementType == Elements.triangle2)
    xref = [0 0;1 0;0 1;0.5 0;0.5 0.5;0 0.5]';
    N = Nodal2D(Operators.curl);
    CB_plot = cell( size(xref,2), 1);
    
    %arrays for plotting
    X_plot = zeros(6, Ne); 
    Y_plot = zeros(6, Ne);
    Physt_plot = zeros(6, Ne);
    Peddy_plot = zeros(6, Ne);
    Pexcess_plot = zeros(6, Ne);
    for k = 1:numel(CB_plot)
        Bx = zeros(Ne, Nsamples);
        By = zeros(Ne, Nsamples);
        for k_shape = 1:N_shape
            X_shape = A(sim.msh.t(k_shape,:), :);
            Nhere = N.eval(k_shape, xref(:,k), msh, 1:Ne);
            
            Bx = Bx + bsxfun(@times, Nhere(1,:)', X_shape );
            By = By + bsxfun(@times, Nhere(2,:)', X_shape );
        end
        %fft
        Cx = fft(Bx, [], 2)/L*2;
        Cy = fft(By, [], 2)/L*2;
        CB_plot{k} = abs(Cx).^2 + abs(Cy).^2;
        
        %plotting points
        X_plot(k,:) = msh.p(1, msh.t(k, :));
        Y_plot(k,:) = msh.p(2, msh.t(k, :));
    end
end

%computing losses
Physt = zeros(Ne, 1);
Peddy = zeros(Ne, 1);
Pexcess = zeros(Ne, 1);
Ptot = 0;


%mats = 13;
mats = unique(msh.matel);
for k = 1:numel(mats)
    [~, coeffs, density] = get_defaultMaterials( mats(k) );
    %coeffs = [1.41 0.061 0];
    if ~any(coeffs)
        continue;
    end
    
    ch = coeffs(1); ce = coeffs(2); cex = coeffs(3);
    els = find(msh.matel == mats(k));
    
    for k_quad = 1:N_quad
        %hysteresis loss density
        temp = bsxfun(@times, CB2{k_quad}(els,1:Nrel), ch*freqRatios(1:Nrel));
        Physt(els) = Physt(els) + Jc.w(k_quad) * sum(temp, 2);
        
        %eddy-current loss density
        temp = bsxfun(@times, CB2{k_quad}(els,1:Nrel), ce*freqRatios(1:Nrel).^2);
        Peddy(els) = Peddy(els) + Jc.w(k_quad) * sum(temp, 2); %loss density
        
        %excess loss density
        temp = bsxfun(@times, CB2{k_quad}(els,1:Nrel).^0.75, cex*freqRatios(1:Nrel).^1.5);
        Pexcess(els) = Pexcess(els) + Jc.w(k_quad) * sum(temp, 2); %loss density
        
    end
    
    Ptot = Ptot + ...
        sum(bsxfun(@times, Physt(els)+Peddy(els)+Pexcess(els), ...
        density*Jc.DETF(els)'))*sim.dims.leff*sim.msh.symmetrySectors;
    
    %values for plotting?
    if pretty_plotting && (msh.elementType == Elements.triangle2)
        for kplot = 1:numel(CB_plot)
            %hysteresis loss density
            temp = bsxfun(@times, CB_plot{kplot}(els,1:Nrel), ch*freqRatios(1:Nrel));
            Physt_plot(kplot, els) = Physt_plot(kplot, els) + sum(temp, 2)';
            
            %eddy-current loss density
            temp = bsxfun(@times, CB_plot{kplot}(els,1:Nrel), ce*freqRatios(1:Nrel).^2);
            Peddy_plot(kplot, els) = Peddy_plot(kplot, els) + sum(temp, 2)'; %loss density
            
            %excess loss density
            temp = bsxfun(@times, CB_plot{kplot}(els,1:Nrel).^0.75, cex*freqRatios(1:Nrel).^1.5);
            Pexcess_plot(kplot, els) = Pexcess_plot(kplot, els) + sum(temp, 2)'; %loss density
        end
    end
end

%Physt_plot
%X_plot

if plotting_on
    if pretty_plotting && (msh.elementType == Elements.triangle2)
        ord = [1 4 2 5 3 6];
        
        figure(8); clf; hold on; box on; axis equal tight;
        patch(X_plot(ord,:), Y_plot(ord,:), Physt_plot(ord,:), 'linestyle', 'none', ...
            'FaceColor','interp'); 
        colormap('jet'); colorbar; caxis([0 [0 1]*caxis']);
        title('Hysteresis loss density (W/kg)')
        
        figure(9); clf; hold on; box on; axis equal tight;
        patch(X_plot(ord,:), Y_plot(ord,:), Peddy_plot(ord,:), 'linestyle', 'none', ...
            'FaceColor','interp'); 
        colormap('jet'); colorbar; caxis([0 [0 1]*caxis']);
        title('Eddy-current loss density (W/kg)')
        
        figure(10); clf; hold on; box on; axis equal tight;
        patch(X_plot(ord,:), Y_plot(ord,:), Pexcess_plot(ord,:), 'linestyle', 'none', ...
            'FaceColor','interp'); 
        colormap('jet'); colorbar; caxis([0 [0 1]*caxis']);
        title('Excess loss density (W/kg)')
        
    else
        figure(8); clf; hold on; box on; axis equal tight;
        msh_fill(msh, [], 2*Physt, 'linestyle', 'none'); colormap('jet');
        colorbar; caxis([0 [0 1]*caxis']);
        title('Hysteresis loss density (W/kg)')

        figure(9); clf; hold on; box on; axis equal tight;
        msh_fill(msh, [], 2*Peddy, 'linestyle', 'none'); colormap('jet');
        colorbar; caxis([0 [0 1]*caxis']);
        title('Eddy-current loss density (W/kg)')

        figure(10); clf; hold on; box on; axis equal tight;
        msh_fill(msh, [], 2*Pexcess, 'linestyle', 'none'); colormap('jet');
        colorbar; caxis([0 [0 1]*caxis']);
        title('Excess loss density (W/kg)')
    end
end

end