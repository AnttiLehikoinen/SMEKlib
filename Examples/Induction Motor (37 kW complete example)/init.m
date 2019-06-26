addpath(genpath('..\..\SMEKlib'));

%setting dimensions
dim = struct();
dim.gmsh_path = 'E:/Software/Work/gmsh43';

%winding and other main dimensions
dim.p = 2; %number of pole pairs
dim.leff = 0.249; %core length
dim.symmetrySectors = 4; %number of sectors in simulation
dim.N_series = 12; %number of turns per coil
dim.l_halfCoil = 0.520; %length of half a coil (for EW resistance)
dim.a = 2; %number of parallel paths
dim.N_layers = 1; %number of winding layers
dim.A_ring = 520.00E-06; %cross-sectional area of end-ring
dim.D_ring = 154.00E-03; %diameter of end-ring
dim.connection_stator = defs.star; %stator connection
dim.Lew = 0.1e-3; %end-winding inductance

dim.type_statorWinding = defs.stranded; %stator winding type
dim.fillingFactor = 0.45; %stator filling factor
dim.type_rotorWinding = defs.cage; %rotor winding type

dim.sigma_stator = 45e6; %stator winding conductivity
dim.sigma_rotor = 35.5E+06 *(230.0 + 20)/(230.0 + 80); %rotor winding conductivity

%stator dimensions
dim.Sout = 310.0e-3/2; %outer radius of stator
dim.Sin = 200.0E-03 /2; %inner radius of stator
dim.Qs = 48;
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
dim.Rout = 198.40E-03 / 2;
dim.Rin = 70.00E-03  / 2;
dim.Qr = 40;
dim.RM = 4;

dim.h_bridge_r = 0.7e-3; %bridge thickness
dim.w_bar_out = 6e-3; %outer bar width
dim.w_bar_mid = 2.5e-3; %mid-bar width

dim.h_bar_out = 16.1e-3; %outer bar-thingy height
dim.h_bar_in = 17.8e-3; %inner bar-thingy height