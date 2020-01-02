
ts = linspace(0, 1/pars.f, 5000);

U1 = sqrt(2)*230; %amplitude of fundamental
Ubus = 400; %DC-bus voltage
fs = 1000; %switching frequency

U = ThreePhasePWM(U1, pars.f, ts, Ubus, fs);

%rms(U,2)
temp = fft(U(1,:));
abs(temp(2))/numel(ts)*2

figure(10); clf; hold on;
plot(ts, U(1,:), 'b');
plot(ts, U1*cos(2*pi*50*ts-0*2*pi/3), 'b');

temp2 = fft(U1*cos(2*pi*50*ts-0*2*pi/3));
abs(temp2(2))/numel(ts)*2
