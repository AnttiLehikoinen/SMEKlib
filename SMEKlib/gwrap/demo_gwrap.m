%a mini-minimal demo for the gwrap wrapper

gmsh_path = 'D:\Large Installations\Gmsh\'; %CHANGE THIS
gwrap_init(gmsh_path);

gwrap_addpolygon([0.25 0.45 0.45 0.25], [0.25 0.25 0.5 0.5], 'Some_rectangle');
gwrap_addcircle(0.75, 0.5, 0.25, 100, 'Some_circle');
gwrap_addpolygon([0 1.05 1.05 0], [0 0 1 1], 'OuterBoundary');

[t, p, t_inEntity, Name2id] = gwrap_finalize();

figure(1); clf; hold on; box on;
triplot(t(:,:)', p(1,:), p(2,:), 'c');
triplot(t(:, t_inEntity==Name2id('Some_rectangle'))', p(1,:), p(2,:), 'b');
triplot(t(:, t_inEntity==Name2id('Some_circle'))', p(1,:), p(2,:), 'r');
legend('Other stuff', 'Some rectangle', 'Some circle');
axis equal;
