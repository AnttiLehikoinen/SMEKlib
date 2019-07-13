%init Run a complete analysis workflow on a fictional 12-phase 37 kW induction motor, 
% from meshing to post-processing. As there is not "official" support for 
% multi-phase machines, some manual setting is needed for the supply
% voltage vectors and stator winding matrices.
%
% By default, the free open-source software gmsh is used for meshing the geometry; 
% you can download yours from
%   http://gmsh.info/#Download
% gmsh is covered by its own separate license. At least until version 4.3, no installation
% is required; it is sufficient to just download and extract the files.
%
% Alternatively, you can generate your own mesh with the tool of your choice and import it.
%
% (c) 2019 Antti Lehikoinen / Smeklab Ltd
%   www.smeklab.com
%   antti@smeklab.com

addpath(genpath('../../SMEKlib'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setting dimensions
dim = struct();
dim.gmsh_path = 'C:/Users/Antti/Documents/gmsh'; %CHANGE THIS ACCORDINGLY

%winding and other main dimensions
dim.p = 2; %number of pole pairs
dim.leff = 0.249; %core length
dim.symmetrySectors = 4; %number of sectors in simulation
dim.N_series = 12; %number of turns per coil
dim.l_halfCoil = 0.520; %length of half a coil (for EW resistance)
dim.a = 1; %number of parallel paths
dim.N_layers = 1; %number of winding layers
dim.A_ring = 520.00E-06; %cross-sectional area of end-ring
dim.D_ring = 154.00E-03; %diameter of end-ring
dim.connection_stator = defs.delta; %stator connection

%setting multi-phase stuff
dim.Qs = 48; %number of stator slots
slot_angles = linspace(0.5, dim.Qs-0.5, dim.Qs)'*2*pi/dim.Qs * dim.p;
slot_angles = slot_angles( 1:(dim.Qs/dim.symmetrySectors) ); %slot angles in el. rads

Uampl = 400/sqrt(3)*sqrt(2) / 8;
Uh = Uampl * exp(1i*-slot_angles); %phase voltage phasors
Ut = @(t)( Uampl * cos( bsxfun(@minus, 2*pi*50*t, slot_angles) ) ); %voltage function

%winding matrix
Ls = kron( [1;-1;1;-1], eye(12) );
dim.Ls = Ls;

dim.type_statorWinding = defs.stranded; %stator winding type
dim.fillingFactor = 0.45; %stator filling factor
dim.type_rotorWinding = defs.cage; %rotor winding type

dim.sigma_stator = 45e6; %stator winding conductivity
dim.sigma_rotor = 35.5E+06 *(230.0 + 20)/(230.0 + 80); %rotor winding conductivity

%stator dimensions
dim.Sout = 310.0e-3/2; %outer radius of stator
dim.Sin = 200.0E-03 /2; %inner radius of stator
dim.delta = 0.8e-3; %airgap
dim.SM = 4; %stator core material

%stator slot dimensions
dim.hslot_s = 23.9e-3; %total slot height, airgap to slot bottom
dim.htt_taper_s = 1e-3; %slot opening, straight part before taper

dim.ws_o = 6.5e-3; %slot width, opening side
dim.ws_b = 8.8e-3; %slot width, bottom side
dim.w_slotOpening_s = 3.5e-3; %slot opening width

dim.htt_s = dim.hslot_s - dim.ws_b/2 - 17.5e-3; %slot opening, total height

%rotor dimensions
dim.Rout = dim.Sin - dim.delta; %outer radius of rotor
dim.Rin = 70.00E-03  / 2; %inner (shaft) radius of rotor
dim.Qr = 40; %number of rotor bars
dim.RM = 4; %rotor core material

dim.h_bridge_r = 0.7e-3; %bridge thickness (outer bar surface to airgap)
dim.w_bar_out = 6e-3; %outer bar width
dim.w_bar_mid = 2.5e-3; %mid-bar width

dim.h_bar_out = 16.1e-3; %outer bar-thingy height (distance from outer bar surface to inner bar top-arc center)
dim.h_bar_in = 17.8e-3; %inner bar-thingy height (distance from inner bar outer-arc center to inner-arc center)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running analysis scripts
mesh_stator; %generate mesh data for stator
mesh_rotor; %generate mesh data for rotor
parse_mesh; %create mesh object

run_simulation_multiphase; %run simulation and compute losses