classdef ThreePhaseController < handle
    properties
        R, L, alpha, Tp, TpInv,
        ki, kp, I, Ra
        Iref, tprev,
        Is, Us, ri
        Icorr, w, fs, Ubus, angle_bias
        
        efilt
    end
    
    methods
        function this = ThreePhaseController(Iref, fs, Ubus, pars, dims)
            
            this.R = 0.4656;
            this.L = 0.0122;
            this.alpha = 500;
            this.w = 2*pi*pars.f;
            this.fs = fs;
            this.Ubus = Ubus;
            
            this.kp = this.alpha*this.R;
            this.ki = this.alpha^2 * this.L;
            this.I = [0; 0];
            this.Ra = this.alpha*this.L - this.R;
            %this.Ra = 0*diag([this.Ra 0*this.Ra]);
            
            shift = +1/6;
            angles = -[0 2/3 4/3 0+shift 2/3+shift 4/3+shift]*pi;
            this.Tp = 1/3*[cos(angles); -sin(angles)];
            this.angle_bias = 0.0654;  
            this.TpInv = 3*this.Tp';
            
            %3-phase configuration
            this.Tp = 2*this.Tp(:,1:3);
            this.angle_bias =  0.0654;            
            this.TpInv = 3/2*this.Tp';
            
            %phi_rotor = pi/(2*dims.p);
            %q = dims.q; c = dims.c;
            %phi_stator = pi/2/dims.p + 2*pi/dims.Qs*(0.5*(0.5 + (q-1)/2)+0.5*(0.5+c+(q-1)/2));
            %phi_stator = phi_stator+pi/dims.p;
            %this.angle_bias = dims.p*(phi_rotor-phi_stator);

            
            
            %this.Iref = [0;20*3/2*sqrt(2)];
            this.Iref = Iref;
            this.tprev = 0;
            this.Icorr = [0;0];
            
            Nsamples = numel(pars.ts);
            this.Is = zeros(2, Nsamples);
            this.Us = zeros(3, Nsamples);
            this.ri = 1;
            
            this.efilt = zeros(2, 3);
        end
        
        function e = lowpass(this, e)
           this.efilt = [e this.efilt(:,1:(end-1))];
           e = mean(this.efilt, 2);
        end
        
        function Uout = U(this, Is, rotorAngle, t)
            %ra = 2*rotorAngle  + 1.1781;
            %ra = 2*rotorAngle  + 1.0472; %3-phase
            ra = 4*rotorAngle + this.angle_bias; %3-phase
            
            rotM = [cos(-ra) -sin(-ra);
                sin(-ra) cos(-ra)];
            Idq = rotM*this.Tp*Is;
            
            dt = t - this.tprev; this.tprev = t;
            
            e_orig = this.Iref - Idq;
            %e_orig = this.lowpass(e_orig);
            e = e_orig + this.alpha*this.Icorr;
            
            %anti-windup for correction
            %e_orig = min(5, abs(e_orig)).*sign(e_orig);
            e_orig = e_orig.*(1 -1*(abs(e_orig)>10));
            this.Icorr = this.Icorr + 5*dt*e_orig;
            %disp(this.Icorr);
            
            %uref = this.kp*e + this.ki*this.I - this.Ra*Idq;
            uref = this.kp*e + this.ki*this.I - [this.R this.w*this.L;
                -this.w*this.L this.R]*Idq;
            
            %getting output voltage          
            Uout_ref = this.TpInv*rotM'*uref;
            
            this.Is(:,this.ri) = Idq;
            Uout = this.Upwm(Uout_ref, t);
            this.Us(:,this.ri) = Uout;
            this.ri = this.ri+1;
            
            %output limitation
            %Umax = 750;
            %v = min(Umax*[1; 1], abs(Uout)).*sign(uref);
            %Uout = Uout * min(1, Umax/norm(Uout));
            
            Udq = rotM*this.Tp*Uout;
            
            %integration with anti-windup
            %this.I = this.I + dt*(e + 1/this.kp*(v - uref));
            this.I = this.I + dt*(e + 1/this.kp*(Udq - uref));
            
            %output limitation
            %Umax = 750;
            %v = min(Umax*[1; 1], abs(Uout)).*sign(uref);
            %Uout = Uout * min(1, Umax/norm(Uout));
            
            %realized voltage space vector
            Udq = rotM*this.Tp*Uout;
            
            %integration with anti-windup
            %this.I = this.I + dt*(e + 1/this.kp*(v - uref));
            this.I = this.I + dt*(e + 1/this.kp*(Udq - uref));
            
            return
            figure(10);
            plot( t, Idq(1), 'b.');
            plot( t, Idq(2), 'r.');
            plot( t, mean(this.Is(:, max(1,this.ri-10):(this.ri-1)),2), 'kx');
            
            drawnow;
            
            figure(11);
            plot(t, Uout(1), 'b.');
            plot(t, Uout(2), 'r.');
            plot(t, Uout(3), 'k.');
            drawnow;
        end
        
        function U = Upwm(this, Ud, t)
           %getting instantaneous PWM output
           %U = Ud;            return;
           
           Ubase = this.Ubus / sqrt(3);
           
           %getting control signals
           Uy = [1 -1 0;0 1 -1;1 1 1] \ [Ud(1:2); 0] / Ubase;
           
           %symmetricing
           delta = (max(Uy) + min(Uy))/2;
           Uy = Uy - delta;

           %overmodulating if necessary
           m = max(Uy);
           if m > 1
               Uy = Uy/m;
           end
           
           %triangle comparison
           Utr = sawtooth(2*pi*this.fs*t);
           Upwm = this.Ubus*( Uy >= Utr );

           Upwm = Uy*Ubase;

           U = [1 -1 0;0 1 -1;-1 0 1]*Upwm;
        end
    end
end