%time-stepping simulation with plots

%%{
f = 500;
w = 2*pi*f;

N_stepsPerPeriod = 100;
N_periods = 2;

ts = linspace(0, N_periods/f, N_stepsPerPeriod*N_periods);
Nsamples = numel(ts);
dt = ts(2) - ts(1);

Aprev = zeros(size(msh.edges,2), 1);

Jvecs = zeros(2, Ne, Nsamples);
Bvecs = zeros(1, Ne, Nsamples);

Q = [S_AA(e_free, e_free)+1/dt*M_AA(e_free,e_free) V_Constr(e_free);
    V_Constr(e_free)' 0];


for ks = 2:Nsamples
    disp(['Computing timestep ' num2str(ks)])
    
    
    FL = [M_AA(e_free,e_free)*Aprev(e_free)/dt; Phi*sin(w*ts(ks))];    
    X = Q \ FL;
    
    Anew = 0*Aprev;    
    Anew(e_free) = X(e_free);
    dA = (Anew - Aprev)/dt;
    
    %computing J and B
    for k = 1:3
        Jvecs(:,:,ks) = Jvecs(:,:,ks) + bsxfun(@times, W.eval(k, xref, msh, []), sigma*dA(abs(msh.t2e(k,:)))' );
        Bvecs(:,:,ks) = Bvecs(:,:,ks) + bsxfun(@times, Wcurl.eval(k, xref, msh, []), Anew(abs(msh.t2e(k,:)))' );
    end
    
    %plotting
    figure(4); clf; hold on; box on; axis equal tight;
    quiver(Xplot(1,:), Xplot(2,:), Jvecs(1,:,ks), Jvecs(2,:,ks)); 
    axis([0 0.0025 0 1e-3]);
    drawnow;
    
    Aprev = Anew;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating animation

axl = [0 0.005 -0.1e-3 1.1e-3];

els = find( Xplot(1,:) < (axl(2) + 1e-3) );

h = figure(4); clf; 
steps = 100:2:200;
filename = 'animation2.gif'; n = 1;
for k = 1:numel(steps)
    ks = steps(k);
    
    
    subplot(2,1,1);box on; axis equal tight;
    quiver(Xplot(1,:), Xplot(2,:), Jvecs(1,:,ks), Jvecs(2,:,ks));
    axis(axl);
    title('Current density')
    
    subplot(2,1,2); box on; axis equal tight;
    fill(Xtri(:,els), Ytri(:,els), Bvecs(:,els,ks), 'linestyle', 'none'); 
    colormap('jet'); colorbar; caxis([-2 2])
    axis(axl);
    title('Flux density');
    xlabel(['Time = ' num2str(1e3*ts(ks)) ' ms'])
    
    drawnow;
    
    %saving
    %saving
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.001);
    elseif k < numel(steps)
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.001);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1.5);
    end
end