function U = PWM_1pulse(x, U1, Ubus)
%PWM_1pulse 1-pulse PWM.
%
% U = PWM_1pulse(t, f, U1, Ubus)
% returns the voltage U at angular position x, PWM-modulated with 1 pulse per 
% half-period, with its fundamental (at f Hz) amplitude equal to U1. The
% instantaneous value U(x) is equal to +- Ubus.
%
% (c) 2018 Antti Lehikoinen / Aalto University.

%Determining the duty ratio for Matlab's square function:
% based on the fact that the fundamental of a square wave is equal to
% sin(a/2), where a = length of pulse in radians [0,pi].

r = min(U1 /(4/pi*Ubus), 1);
a = asin(r);
d = a/(pi/2)*100;

auxfun = @(x)( square(x + pi*d/200, d/2) );

U = (0.5+0.5*auxfun(x) + ...
    -0.5-0.5*auxfun(x-pi))*Ubus;

end