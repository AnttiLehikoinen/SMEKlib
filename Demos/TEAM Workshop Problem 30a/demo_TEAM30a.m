%simulates the TEAM workshop problem 30a
% 
% (c) 2017 Antti Lehikoinen / Aalto University and Jonathan Velasco 

gmsh_path = 'C:\Antti\Software\gmsh'; %CHANGE THIS

addpath(genpath('..\..\SMEKlib')); %adding library to the search path

%parameters
mur_stator = 30; %relative permeability
mur_rotor = 30;
sigma_rotorSteel = 1.6e6; %conductivity
sigma_aluminum = 3.72e7;
I = 2045.175; %stator current

fs = 60;
wms = linspace(0, 1200, 49); %range of rotor speeds analysed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining geometry

%meshing and extracting geometry
geo_path = 'Problem30a.geo';

gwrap_mesh(gmsh_path, geo_path);
[t, p, t_inEntity, Name2id] = gwrap_loadmesh(geo_path);
disp(Name2id.keys)

% parsing data:
%fixing stator conductor order for counter-clockwise rotation
statorConductors = cell(1,6);
statorConductors{1} = find(t_inEntity == Name2id('A'));
statorConductors{6} = find(t_inEntity == Name2id('-C'));
statorConductors{5} = find(t_inEntity == Name2id('B'));
statorConductors{4} = find(t_inEntity == Name2id('-A'));
statorConductors{3} = find(t_inEntity == Name2id('C'));
statorConductors{2} = find(t_inEntity == Name2id('-B'));

%other elements in non-air domains
statorIron = find(t_inEntity == Name2id('Stator Steel'));
rotorIron = find(t_inEntity == Name2id('Rotor Steel'));
rotorCoat = find(t_inEntity == Name2id('Rotor Aluminum'));

%winding configuration
W = windingConfiguration_1(1, 1); %winding configuration matrix
L = statorConnectionMatrix(W, 1, 1);

%air-gap elements
el_ag = [find(t_inEntity == Name2id('Air Gap to Rotor')) find(t_inEntity == Name2id('Air Gap to Stator'))];

%initializing mesh struct
msh = inittri(p, t);
msh.rotel = [rotorIron rotorCoat];
msh.bandData = initializeBandData(msh, [], rotorCoat, msh.t(:, el_ag));

%number of nodes & elements, plus boundary nodes
Np = size(msh.p, 2);
Ne = size(msh.t, 2);
n_Dirichlet = [unique( msh.edges(:, ~msh.e2t(2,:)) )' find((abs(msh.p(1,:))<1e-4) & (abs(msh.p(2,:))<1e-4))];
n_free = setdiff(1:Np, n_Dirichlet);

%plotting
figure(1); clf; hold on; box on;
msh_fill(msh, statorIron, [1 1 1]*0.4);
msh_fill(msh, rotorIron, [1 1 1]*0.4);
msh_fill(msh, rotorCoat, 'c');
msh_fill(msh, cell2mat(statorConductors), 'r');
msh_fill(msh, find(t_inEntity == Name2id('Air Between Coils')), 'y');
msh_plot(msh, n_Dirichlet, 'ko');

%plotting the 6th stator coil just because
msh_fill(msh, statorConductors{6}, 'g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling matrices
tic;

%setting reluctivities per element
mu0 = pi*4e-7;
nu = 1/mu0 * ones(1,Ne);
nu(statorIron) = 1/(mur_stator*mu0);
nu(rotorIron) = 1/(mur_rotor*mu0);

%assembling stiffness matrix
S_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', nu, [], msh, []);
S = sparseFinalize(S_struct, Np, Np);

%assembling mass matrix in two parts
M_struct = assemble_matrix('', 'nodal', '', 'nodal', sigma_rotorSteel, rotorIron, msh, []);
M_struct = assemble_matrix('', 'nodal', '', 'nodal', sigma_aluminum, rotorCoat, msh, M_struct);
M = sparseFinalize(M_struct, Np, Np);

%supply matrix
JF_struct = [];
for k = 1:6
    JF_struct = assemble_vector('', 'nodal', 1, k, statorConductors{k}, msh, JF_struct);
end
JC = sparseFinalize(JF_struct, Np, 6);
cA = sum(JC*eye(6),1); %conductor areas
JC = bsxfun(@times, JC, 1./cA);

%setting load vector
Iphase = [1; exp(-1i*2*pi/3); exp(-1i*4*pi/3)]*I;
FL = JC*L*Iphase;

disp(['Matrices assembled in ' num2str(toc) ' seconds.']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving
tic;
As = zeros(Np, numel(wms));
slips = (2*pi*fs - wms)/(2*pi*fs);
for kw = 1:numel(wms)
    slip = slips(kw);
    
    Q = S + 1i*2*pi*fs*slip*M;
    As(n_free, kw) = Q(n_free, n_free) \ FL(n_free);
end
disp([num2str(numel(wms)) ' operating points solved in ' num2str(toc) ' seconds.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing

%induced back-emfs
Es_coil = 1i*2*pi*fs*JC'*As;

figure(2); clf; hold on; box on;
hE = plot(wms, abs(Es_coil), 'b');
title('Induced voltages per turn and coil-side')
ylabel('Voltage (V)')
xlabel('Speed (rad/s)');

%computing current densities in rotor
J_coat = zeros(Np, numel(wms)); J_iron = zeros(Np, numel(wms));
J_coat( msh.t(:,unique(rotorCoat)),: ) = -1i*2*pi*fs*sigma_aluminum*bsxfun(@times, As(msh.t(:,unique(rotorCoat)),:), slips);
J_iron( msh.t(:,unique(rotorIron)),: ) = -1i*2*pi*fs*sigma_rotorSteel*bsxfun(@times, As(msh.t(:,unique(rotorIron)),:), slips);

%mass-type matrices for integrating the loss density in rotor
M_coat = sparseFinalize(assemble_matrix('', 'nodal', '', 'nodal', 1, rotorCoat, msh, []), ...
    Np, Np);
M_iron = sparseFinalize(assemble_matrix('', 'nodal', '', 'nodal', 1, rotorIron, msh, []), ...
    Np, Np);

P_coat = zeros(1, numel(wms)); P_iron = P_coat;
for kw = 1:numel(wms)
    P_coat(kw) = real( (J_coat(:,kw))'*M_coat*(J_coat(:,kw)/sigma_aluminum) );
    P_iron(kw) = real( J_iron(:,kw)'*M_iron*(J_iron(:,kw)/sigma_rotorSteel) );
end

%plotting current ant flux densities at a specified operating point
kw = 9; %200 rad/s

figure(3); clf; hold on; box on;
drawCurrentDensity(msh, real(J_coat(:,kw)), rotorCoat );
drawCurrentDensity(msh, real(J_iron(:,kw)), rotorIron ); colorbar;
title(['Rotor current density at ' num2str(wms(kw)) ' rad/s']);

figure(4); clf; hold on; box on;
subplot(2, 1, 1); hold on; box on;
plot(wms, P_coat, 'b');
title('Losses in aluminum coat (W)');

subplot(2, 1, 2); hold on; box on;
plot(wms, P_iron, 'r');
title('Losses in rotor iron (W)');
xlabel('Speed (rad/s)');

%plot Comsol results if available
try
    Data_Comsol = csvread('Data.csv', 5, 0);
    
    figure(2);
    hE_Comsol = plot(Data_Comsol(:,1), abs(Data_Comsol(:, 2:7)), 'rv');
    legend([hE(1) hE_Comsol(2)], 'Matlab', 'Comsol');
    
    figure(4);
    subplot(2,1,1);
    plot(Data_Comsol(:,1), Data_Comsol(:,8), 'bv');
    legend('Matlab', 'Comsol', 'Location', 'Best');
    subplot(2,1,2);
    plot(Data_Comsol(:,1), Data_Comsol(:,9), 'rv');
    legend('Matlab', 'Comsol', 'Location', 'Best');
catch
    print('No Comsol simulation results found!');
end

figure(5); clf; hold on; box on;
drawFluxDensity(msh, real(As(:,kw)), 'LineStyle', 'none'); colormap('jet'); 
colorbar; caxis([0 0.22]);
drawFluxLines(msh, real(As(:,kw)), 20, 'k');
title(['Flux density and flux lines at ' num2str(wms(kw)) ' rad/s']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performing time-stepping analysis because WE CAN
disp('Performing time-stepping analysis...');

Nperiods = 2;
N_stepsPerPeriod = 400;

dt = 1/fs / N_stepsPerPeriod; %time-step length

%re-assembling stiffness matrix for non-airgap elements
Sstruct_td = assemble_matrix('grad', 'nodal', 'grad', 'nodal', nu, setdiff(1:Ne, el_ag), msh, []);
S = sparseFinalize(Sstruct_td, Np, Np);

%allocating arrays
Nsamples = ceil( Nperiods/fs / dt );
tsamples = (0:(Nsamples-1))*dt;
Asamples = zeros(Np, Nsamples); Asamples(:,1) = real(As(:,kw))*sqrt(2);
Esamples = zeros(3, Nsamples);

%performing time-stepping with backward-Euler
tic;
for ksample = 2:Nsamples
    theta = 2*pi*fs*tsamples(ksample);
    FL = JC*L*[cos(theta); cos(theta - 2*pi/3);  cos(theta-4*pi/3)]*I*sqrt(2) + ...
        1/dt*M*Asamples(:,ksample-1);
    
    Q = S + 1/dt*M + get_MovingBandMatrix(wms(kw)*tsamples(ksample), msh);
    
    Asamples(n_free, ksample) = Q(n_free,n_free) \ FL(n_free);
    Esamples(:, ksample) = L'*JC'*(Asamples(:,ksample)-Asamples(:,ksample-1))/dt;
end
disp([num2str(Nsamples) ' time-steps simulated in ' num2str(toc) ' seconds.']);

figure(6); clf;
plot(tsamples, Esamples');

%plotting current density in the rotor
kt = 800; %time-step to plot at

J_coat_ts = zeros(Np, 1); J_iron_ts = zeros(Np, 1);
J_coat_ts( msh.t(:,unique(rotorCoat)) ) = -sigma_aluminum*...
    (Asamples(msh.t(:,unique(rotorCoat)),kt) - Asamples(msh.t(:,unique(rotorCoat)),kt-1))/dt;
    
J_iron_ts( msh.t(:,unique(rotorIron)) ) = sigma_rotorSteel*...
    (Asamples(msh.t(:,unique(rotorIron)),kt) - Asamples(msh.t(:,unique(rotorIron)),kt-1))/dt;

figure(7); clf; hold on; box on;
drawCurrentDensity(msh, J_coat_ts, rotorCoat );
drawCurrentDensity(msh, J_iron_ts, rotorIron ); colorbar;axis equal tight;