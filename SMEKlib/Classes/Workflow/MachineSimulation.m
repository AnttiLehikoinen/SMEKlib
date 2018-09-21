classdef MachineSimulation < handle
    %MachineSimulation Base class for analysis of rotating machines.
    % 
    % (c) 2017 Antti Lehikoinen / Aalto University
    
    properties
        msh, dims, results, Ne, Np, matrices, nu_fun, nu_struct, misc
    end
    
    methods
        function this = MachineSimulation(msh, dims)
            %constructor
            this.msh = msh;
            this.Np = size(msh.p, 2);
            this.Ne = size(msh.t, 2);
            this.dims = dims;
            this.misc = struct();
            
            this.matrices = struct('P', [], 'W', [], 'Ls', [], 'Lr', [], ...
                'Ms', [], 'Mr', [], ...
                'Cs', [], 'Cr', [], ...
                'DRs', [], 'DRr', [], 'Zew_s', [], 'Zew_r', []);
            this.results = struct('Xh', [], 'Xt', []);
            
            this.nu_struct = initialize_reluctivityStruct_interp1(msh, true);
            this.nu_fun = @(B)( calculate_reluctivity(B, this.nu_struct) );
            
            this.setBoundaryMatrix();
            this.setStatorCircuitMatrices();
            this.setRotorCircuitMatrices();
            this.setLoadVector();
        end
        
        function this = setBoundaryMatrix(this, varargin)
            if numel(varargin)
                %TODO
            else
                n_dir = this.msh.namedNodes.get('Dirichlet');
                np_master = this.msh.namedNodes.get('Periodic_master');
                np_slave = this.msh.namedNodes.get('Periodic_slave');
                P_data = {[np_slave; np_master; this.msh.periodicityCoeff*ones(1, numel(np_master))], ...
                [n_dir; zeros(2, numel(n_dir))]};
                this.matrices.P = assemble_TotalMasterSlaveMatrix(this.Np, P_data, []);
            end
        end
        
        function this = setLoadVector(this)
            %setting PM sources, if any
            PMs = this.msh.namedElements.get('PMs');
            if isempty(PMs)
                this.matrices.F = sparse(this.Np, 1);
                return
            end
            Fc = MatrixConstructor();
            for k = 1:size(PMs,2)
                Fc.assemble_vector(Nodal2D(Operators.curl), 1, PMs{1,k}, PMs{2,k}, this.msh);
            end
            this.matrices.F = Fc.finalize(this.Np, 1);
        end
        
        function this = setStatorCircuitMatrices(this)
            this = sim_setStatorCircuitMatrices(this);
        end
        
        function this = setRotorCircuitMatrices(this, varargin)
            this = sim_setRotorCircuitMatrices(this);
        end
        
        function this = run_harmonic(this, varargin)
            %Time-harmonic analysis.
            %
            % Call syntax:
            %   this.run_harmonic(pars)
            %   this.run_harmonic(pars, nu_fun)
            %   this.run_harmonic(pars, BH_fun)
            %
            % Method runs nonlinear "time-harmonic analysis", i.e. solves
            % for the sine and cosine terms of the fundamental-only Fourier
            % series. Thus, the cosine term corresponds to the real part of
            % the corresponding phasor, while the sine is the NEGATIVE
            % imaginary part. For this reason, complex voltages in pars.U
            % must be conjugated to obtain correct results.
            
            this = sim_runTimeHarmonicSimulation(this, varargin{:});
        end               
        function this = init(this, varargin)
            this = sim_initialConditions_CN(this, varargin{:});
        end
        function this = run_timestepping(this, varargin)
            %this = sim_runTimeSteppingSimulation(this, varargin{:});
            this = sim_runTimeSteppingSimulation_CN(this, varargin{:});
        end
        function this = run_static(this, varargin)
            this = sim_runStaticSimulation(this, varargin{:});
        end
        
        function [] = fluxplot(this, step, pars, varargin)
            %Flux density plot.
            % 
            % Call syntax:
            %   [] = this.fluxplot(step, pars)
            %       for time-step 'step'.
            %
            %   [] = this.fluxplot(-1, pars)
            %       for harmonic analysis results.
            
            if ~isempty(pars.slip)
                slip = pars.slip;
            else
                slip = this.dims.slip;
            end
            
            if step == -1
                A = this.results.Xh(1:this.Np, 1);
                step = 1;
                rotorangle = (step-1)*(1-slip)*(2*pi*pars.f/this.dims.p) * ...
                    (1/pars.f) / pars.N_stepsPerPeriod;
            elseif numel(varargin)>0 && strcmp(varargin{1}, 'static')
                A = this.results.Xs(1:this.Np, step);
                if isempty(pars.rotorAngle)
                    rotorangle = 0;
                else
                    rotorangle = pars.rotorAngle;
                    rotorangle = rotorangle(step);
                end
            else
                A = this.results.Xt(1:this.Np, step);
                rotorangle = (step-1)*(1-slip)*(2*pi*pars.f/this.dims.p) * ...
                    (1/pars.f) / pars.N_stepsPerPeriod;
            end
            
            drawFluxDensity(this.msh, A, rotorangle, 'LineStyle', 'none'); 
            colormap('jet'); colorbar; caxis([0 2]);
            drawFluxLines(this.msh, A, 16, rotorangle, 'k');
            axis(this.dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);
            title('Flux lines and flux density (T)');
        end
        
        function I = Is(this)
            indI = this.Np + this.results.Nu_s + this.results.Nu_r + ...
                (1:this.results.Ni_s);
            if this.dims.connection_stator == defs.star
                I = [1 0;0 1;-1 -1] * this.results.Xt(indI,:);
            else
                I = this.results.Xt(indI,:);
            end
        end
        function I = Ish(this)
            %Ish current from time-harmonic analysis
            indI = this.Np + this.results.Nu_s + this.results.Nu_r + ...
                (1:this.results.Ni_s);
            Nvars = size(this.results.Xh, 1)/2;
            I = this.results.Xh(indI) + 1i*this.results.Xh(Nvars+indI);
            if this.dims.connection_stator == defs.star
                I = [1 0;0 1;-1 -1] * I;
            end
        end
            
    end
end