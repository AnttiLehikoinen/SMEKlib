function h = get_periodicityFactor(N_sec, p)
%get_periodicityFactor returns the periodicity factor.
% 
% Call syntax:
% h = get_periodicityFactor(N_sec, p)
% 
% Use real(h) for time-stepping analysis.
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

h = exp(1i * p * 2*pi/N_sec);

if imag(h) < 1e-4
    h = real(h);
end

end