%Performs a nonlinear time-stepping simulation
% 
% Run setup.m first to load mesh etc.
%
% (c) 2018 Antti Lehikoinen / Aalto University

%%{
f = 500;
w = 2*pi*f;

%simulation parameters
N_stepsPerPeriod = 100;
N_periods = 2;

%timestamps
ts = linspace(0, N_periods/f, N_stepsPerPeriod*N_periods);
Nsamples = numel(ts);
dt = ts(2) - ts(1);

%arrays for unknowns and results
Aprev = zeros(size(msh.edges,2), 1);
Jvecs = zeros(2, Ne, Nsamples);
Bvecs = zeros(1, Ne, Nsamples);

%material function
msh.matel = 2*ones(1, size(msh.t,2));
nu_struct = initialize_reluctivityStruct_interp1(msh, true);
nu_fun = @(B)( calculate_reluctivity(B, nu_struct) );
dH_fun = @(B)( matfun_1D(B, nu_fun) );

%constant part of Jacobian
Q = [1/dt*M_AA(e_free,e_free) V_Constr(e_free);
    V_Constr(e_free)' 0];

Jc = JacobianConstructor(msh, Wcurl, Wcurl, true); %constructor for Jacobian
X = zeros(size(msh.edges,2)+1, 1);

figure(10); clf; hold on; box on;
BH = nu_struct.Bnu{1}; B = sqrt(BH(:,1));
nu = BH(B<1.5,2);
B = B( B<1.5 );
plot(B.*nu, B, 'k');
xlabel('H (A/m)');
ylabel('B (T)');
for ks = 2:Nsamples
    disp(['Computing timestep ' num2str(ks)]);   
    
    FL = [M_AA(e_free,e_free)*Aprev(e_free)/dt; Phi*sin(w*ts(ks))];
    
    for kiter = 1:15
        [J, res] = Jc.eval(X, dH_fun );
        
        res_tot = res + Q*X - FL;
        
        resnorm = norm(res_tot);
        disp([' Newton iteration ' num2str(kiter) ', residual norm ' num2str(resnorm)]);        
        if resnorm < 1e-6
            disp('   Converged.');
            break;
        end
        
        X = X - (Q+J) \ res_tot;
    end
    
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
    title(['B = ' num2str( Phi*sin(w*ts(ks))/(1e-3*15e-3) ) ' T']);
    drawnow;
    
    figure(10);
    plot(-X(end), FL(end)/(1e-3*15e-3), 'bo');
    %plot(FL(end)/(1e-3*15e-3), 1./Anew(end), 'bo');
    
    drawnow;
    pause(0.05);
    
    Aprev = Anew;
end
%return
%}

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating animation
filename = 'animation_nonlinear.gif';

axl = [0 0.005 -0.1e-3 1.1e-3];

els = find( Xplot(1,:) < (axl(2) + 1e-3) );

h = figure(4); clf; 
steps = 200:2:400;
for k = 1:numel(steps)
    ks = steps(k);
    
    
    subplot(2,1,1);box on; axis equal tight;
    quiver(Xplot(1,:), Xplot(2,:), Jvecs(1,:,ks), Jvecs(2,:,ks));
    axis(axl);
    title('Current density')
    
    subplot(2,1,2); box on; axis equal tight;
    fill(Xtri(:,els), Ytri(:,els), Bvecs(:,els,ks), 'linestyle', 'none'); 
    colormap('jet'); colorbar; caxis([-3 3])
    axis(axl);
    %title('Flux density');
    title(['Flux density (' num2str( Phi*sin(w*ts(ks))/(1e-3*15e-3) ) ' T average)']);
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