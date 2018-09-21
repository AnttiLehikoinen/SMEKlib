classdef SimulationParameters < dynamicprops
    %simulation parameters container
    properties
        U, f, N_stepsPerPeriod, N_periods, slip, misc, rotorAngle,
        Is, rotorDisplacement,
        alpha2, phi0
    end
    
    methods
        function this = SimulationParameters(varargin)
            %default parameters
            this.U = 400;
            this.f = 50;
            this.N_stepsPerPeriod = 200;
            this.N_periods = 2;
            this.slip = [];
            this.rotorAngle = 0;
            this.Is = [];
            this.rotorDisplacement;
            this.phi0 = 0; %voltage phase angle
            this.alpha2 = 1.1; %weight factor for k+1 step. 
                %2=implicit Euler; 1 = Crank-Nicolson
            
            this.misc = struct('Info', 'Miscellaneous parameters etc.');
            
            if mod(numel(varargin), 2)
                error('Invalid number of input arguments.');
            end
            
            for k = 1:2:numel(varargin)
                try
                    this.(varargin{k}) = varargin{k+1};
                catch
                    %non-default properties
                    this.misc.(varargin{k}) = varargin{k+1};
                end
            end
        end
        
        function p = miscpar(this, name)
            if isfield(this.misc, name)
                p = this.misc.name;
            else
                p = [];
            end
        end
        
        function t = ts(this)
            dt = (1/this.f) / this.N_stepsPerPeriod;
            t = dt * (0:(this.N_stepsPerPeriod*this.N_periods-1));
        end
    end
end
                    