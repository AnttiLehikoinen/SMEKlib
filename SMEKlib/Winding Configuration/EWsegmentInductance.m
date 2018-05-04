function L = EWsegmentInductance(r, A, dims)
%EWsegmentInductance inductance of an end-winding segment.
% 
% L = EWsegmentInductance(r, dims), where
%   r = mean radial coordinate of end-winding
%   A = cross-sectional area of the end-winding
%   dims = machine dimensions
%
% Calculation follows the formula given by Lipo, T.A. for a random-wound
% winding.
%
% (c) 2018 Antti Lehikoinen / Aalto University

if ~isfield(dims, 'l_halfCoil')
    L = 0;
    return
end

lew = dims.l_halfCoil - dims.leff;
taup = r *pi/dims.p; %length of circumferential EW component
l_n = (lew-taup)/2; %length of normal component
l_t = 2*pi/dims.Qs * r; %length of circumferential slot-pitch part

%area enclosed
dR = sqrt(A/pi); %equivalent radius (equal area)

%inductance of tangential part
L = pi*4e-7 / pi *(l_t/2) * log(l_n / (2*dR));

end