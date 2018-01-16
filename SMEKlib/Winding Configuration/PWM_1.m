%function U = threePhasePWM(U, f, Ubus, 

f = 50;
fs = 2000;

t = linspace(0, 2/f, 400);
U = 1*sin(2*pi*f*t);
Ubus = 1;

N_oversampling = 1;

%t_samples = linspace(min(t), max(t), numel(t)*N_oversampling);
t_samples = t;

%Uint = interp1(t, U, t_samples);
Uint = U;
Upwm = (-1 + 2*(Ubus*sawtooth(2*pi*fs*t_samples, 0.5) <= Uint )) * Ubus;

%computing correction
dt = t(2)-t(1);
for kt = 2:numel(t)
    U1 = Uint(kt-1) + (Uint(kt) - Uint(kt-1))*dt;
    U2 = Upwm(kt-1) + (Upwm(kt) - Upwm(kt-1))*dt;
    
    %Upwm(kt) = Upwm(kt) - (U2-U1);
end

figure(10); clf; hold on;

%plot(t, U, 'b');

%plot(t_samples, Ucorr, 'g');

plot(t_samples, Upwm, 'r:')
plot(t_samples, cumsum(Upwm)*dt*2*pi*f/2, 'go-')
plot(t, cumsum(U)*dt*2*pi*f/2, 'kx-')

%axis([0 t(10) -1 2])
