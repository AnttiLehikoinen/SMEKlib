addpath(genpath('..\..\SMEKlib'));

%setting dimensions
dim = struct();
dim.gmsh_path = 'E:/Software/Work/gmsh43';
dim.p = 2;
dim.symmetrySectors = 4;

dim.type_statorWinding = defs.stranded;
dim.type_rotorWinding = defs.cage;
dim.fillingFactor = 0.4;

dim.sigma_stator = 45e6;
dim.sigma_rotor = 35.5E+06 *(230.0 + 20)/(230.0 + 80);

%stator dimensions
dim.leff = 0.249;
dim.Sout = 310.0e-3/2; %outer radius of stator
dim.Sin = 200.0E-03 /2; %inner radius of stator
dim.Qs = 48;
dim.delta = 0.8e-3;
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