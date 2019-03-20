%ag flux density
[tag, pag] = mshc.bandData.t_ag(0);
msh_ag = SimpleMesh(pag, tag);


%finding nodes to plot
angles = linspace(0, pi/4, 1000);
r = 0.55*dim.Sin + 0.45*dim.Rout;
%r = 0.45*dim.Sin + 0.55*dim.Rout;
x = r*cos(angles);
y = r*sin(angles);
[~, els] = pointSources2Elements_faster(x,y,msh_ag); angles = angles(els>0);els = els(els>0);

e_rad = [cos(angles); sin(angles)];
e_circ = [sin(angles); cos(angles)];

rotorAngles = wm*pars.ts;
Nsamples = size(sim.results.Xt, 2);
Babs_all = zeros(numel(angles), Nsamples);
Brad_all = Babs_all;
Bcirc_all = Babs_all;

for k = 1:Nsamples
    
    %computing B in the airgap
    A_ag = mshc.bandData.tag_solution(sim.results.Xt(:,k));
    %A_ag = sim.results.Xt(:,10);
    
    [tag_local, p] = mshc.bandData.t_ag( rotorAngles(k) );
    msh_ag.p = p;
    msh_ag.t = tag_local;
    
    [Babs, Bvec] = calculate_B(A_ag, msh_ag);

    Brad = dotProduct(Bvec(:, els), e_rad);
    Bcirc = dotProduct(Bvec(:, els), e_circ);
    
    %Brad = interp1(angles, Brad(els), mod(angles+rotorAngles(k), pi/4), 'nearest');

    %{
    figure(10); clf; hold on;
    plot(angles/pi*180, Brad(els), 'b' )
    plot(angles/pi*180, Bcirc(els), 'r' )
    plot(angles/pi*180, Babs(els), 'k' )
    title('Flux density')
    drawnow;
    %}
    
    Babs_all(:,k) = Babs(els)';
    Brad_all(:,k) = Brad(els)';
    Bcirc_all(:,k) = Bcirc(els)';    
end



figure(6); clf; hold on; box on;
%drawFluxDensity(mshc, simc.results.Xh, 'LineStyle', 'none'); colormap('jet'); colorbar;% caxis([0 2])
drawFluxDensity(msh_ag, A_ag, 'LineStyle', 'none'); colorbar;
msh_triplot(mshc, find(mshc.matel==2), 'k');

plot(x, y, 'k.-')

figure(7); clf;
surf(Babs_all)