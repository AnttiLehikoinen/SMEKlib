classdef SimulationParameters
    %simulation parameters container
    properties
        U, f, N_stepsPerPeriod, N_periods, slip, misc
    end
    
    methods
        function this = SimulationParameters(varargin)
            %default parameters
            this.U = 400;
            this.f = 50;
            this.N_stepsPerPeriod = 200;
            this.N_periods = 2;
            this.slip = [];
            this.misc = struct('Info', 'Miscellaneous parameters etc.');
            
            if mod(numel(varargin), 2)
                error('Invalid number of input arguments.');
            end
            
            for k = 1:2:numel(varargin)
                try
                    this.(varargin{k}) = varargin{k+1};
                catch
                    error(['Invalid parameter "' num2str(varargin{k}) '"']);
                end
            end
        end
        
        function t = ts(this)
            dt = (1/this.f) / this.N_stepsPerPeriod;
            t = dt * (0:(this.N_stepsPerPeriod*this.N_periods-1));
        end
    end
end
                    