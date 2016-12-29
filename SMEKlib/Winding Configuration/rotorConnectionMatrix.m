function [Lr, Z_ew, L_ew] = rotorConnectionMatrix(Qr, N_sec, p, varargin)
%rotorConnectionMatrix returns the loop matrices for the rotor cage.
% 
% [Lr, Z_ew, L_ew] = rotorCircuitMatrix(Qr, N_sec, p)
% returns the rotor loop matrix Lr, the end-ring impedance matrix Z_ew, and
% the end-ring loop matrix L_ew.
% Use rotorCircuitMatrix(Qr, N_sec, p, Zbe) or
% rotorCircuitMatrix(Qr, N_sec, p, Zbe, Zring) to supply the bar-end and/or
% end-ring segment impedances. Values not given are assumed zero.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%getting/setting bar-end and end-ring segment impedances
if numel(varargin) >= 1
    zbe = varargin{1};
    if numel(varargin) == 2
        zsc = varargin{2};
    else
        zsc = 0;
    end
else
    zbe = 0;
    zsc = 0;
end

%periodicity factor
h = get_periodicityFactor(N_sec, p);

Qr_sec = Qr / N_sec;
Qe_sec = 2*Qr_sec;

%bar-end and end-ring impedance matrices
Zbe = internal_diag(zbe, Qr_sec);
Zsc = internal_diag(zsc, Qe_sec);

Lr = zeros(Qr_sec); 
L_ew = zeros(Qr_sec, Qe_sec);
Lr(1, [1 Qr_sec]) = [1 -h]; 
L_ew(1, [1 Qr_sec+1]) = [-1 1];
for kbar = 2:Qr_sec
    Lr(kbar, [kbar-1 kbar]) = [-1 1];
    L_ew(kbar, [kbar kbar+Qr_sec]) = [1 1];
end

Lr = transpose(Lr);
L_ew = transpose(L_ew);
    
Z_ew = Lr'*Zbe*Lr + L_ew'*Zsc*L_ew;

end

function D = internal_diag(z, Qr)

if numel(z) > 1
    D = diag(z(1:Qr));
else
    D = z * eye(Qr);
end

end