function U = ThreePhasePWM(U1, f, t, Ubus, fs)


theta_vec = bsxfun(@plus, 2*pi*f*t, -[0;2*pi/3;4*pi/3] - pi/6);
theta_s_vec = repmat(2*pi*fs*t, 3, 1);
%theta_s_vec = bsxfun(@plus, 2*pi*fs*t, pi/6*[1;1;1]);

Uaux = ((Ubus*sawtooth(theta_s_vec, 0.5)) <= (2*U1/sqrt(3)*cos(theta_vec)) ) * Ubus;

%figure(11); clf; hold on;
%plot(t, (U1/sqrt(3)*cos(theta_vec(1,:))), 'b');
%plot(t, Ubus*sawtooth(theta_s_vec(1,:),0.5), 'r');
%plot(t, Uaux(1,:), 'r')

U = [1 -1 0;0 1 -1;-1 0 1]*Uaux;

end


