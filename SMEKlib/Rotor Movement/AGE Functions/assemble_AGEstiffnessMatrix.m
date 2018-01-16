function S_AGE = assemble_AGEstiffnessMatrix(N_pp, dims)

ksi = dims.D_si / dims.D_ro; %air-gap radius ratio

ro = dims.D_ro / 2; %outer radius of rotor

%computing the C and D coefficients (Cr^n + Dr^-n) for stator and rotor
%shape functions (radius-dependent part)
CD_s = zeros(N_pp, 2); CD_r = zeros(N_pp, 2);
for k_pp = 1:N_pp
    temp = [1 1;ksi^k_pp ksi^(-k_pp)] \ eye(2);
    CD_r(k_pp,:) = [temp(1,1) temp(2,1)];
    CD_s(k_pp,:) = [temp(1,2) temp(2,2)];
end

%lazy implementation for now
S_AGE = zeros(4*N_pp);
for k = 1:(2*N_pp)
    
    n = floor( (k-1)/2 ) + 1;
    
    % stator-stator terms
    Ci = CD_s(n, 1); Di = CD_s(n,2);  Cj = Ci; Dj = Di;    
    S_AGE(k, k) = (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1))  + ...
        (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        - n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1));
    
    % stator-rotor terms
    Cj = CD_r(n, 1); Dj = CD_r(n, 2);
    S_AGE(k, 2*N_pp + k) = (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1))  + ...
        (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        - n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1));
    S_AGE(2*N_pp + k, k) = S_AGE(k, 2*N_pp + k);
    
    % rotor-rotor terms
    Ci = CD_r(n, 1); Di = CD_r(n, 2);
    S_AGE(2*N_pp + k, 2*N_pp + k) = (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1))  + ...
        (n/2 * Ci*Cj * (ksi^(2*n) - 1) + ...
        - n^2 * (Ci*Dj + Cj*Di) * log(ksi) + ...
        -n/2 * Di*Dj * (ksi^(-2*n) - 1));
end
S_AGE = pi * sparse(S_AGE) / (pi*4e-7);

end