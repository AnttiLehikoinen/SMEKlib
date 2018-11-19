classdef MachineMesh < MeshBase
    %MachineMesh basic machine mesh to be used with the moving band
    %functionality of SMEKlib.
    %
    % Assumes a typical electrical machine, with a circular stator and
    % rotor, and symmetry sector boundaries passing through the origin.
    %
    % Properties in addition to those of MeshBase:
    %   symmetrySectors = number of symmetry sectors
    %   periodicityCoeff = periodicity coefficient
    %   bandData = band data for moving band generation
    %
    % (c) 2017-2018 Antti Lehikoinen / Aalto University, Smeklab Ltd
    
    properties
        symmetrySectors
        periodicityCoeff
        bandData
    end
    
    methods
        function msh = MachineMesh(varargin)
            msh = msh@MeshBase(varargin{:});
            switch size(msh.t,1)
                case 3
                    msh.elementType = Elements.triangle;
                case 6
                    msh.elementType = Elements.triangle2;
                otherwise
                    error('Elementtype not implemented.');
            end
        end
        
        function msh = setSymmetrySectors(msh, N_sec, dims)
            %setSymmetrySectors sets the number of symmetry sectors and
            %the periodicity coefficient.
            %
            % Call syntax setSymmetrySectors(N_of_sectors, dims),
            % where dims.p is the number of pole pairs.
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            
            msh.symmetrySectors = N_sec;
            msh.periodicityCoeff = get_periodicityFactor(N_sec, dims.p);
        end
        
        function rotel = rotel(msh)
            %compatibility hack for the many other SMEKlib functions
            %expecting msh to be a struct with a field named rotel.
            rotel = msh.namedElements.get('rotel');
        end
        function statel = statel(msh)
            statel = msh.namedElements.get('statel');
        end
        
        function msh = setMovingBand(msh, t_ag)
            %setMovingBand initializes moving band data.
            %
            % Call syntax setMovingBand(t_ag) where t_ag is an air-gap
            % triangulation, with indexing referring to the nodes in p.
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            msh.bandData = initializeBandData(msh, msh.statel, msh.rotel, t_ag);
            if msh.elementType == Elements.triangle2 || msh.elementType == Elements.triangle2I
                msh.bandData.shiftTol = 2*msh.bandData.shiftTol;
                msh.bandData.msh_ag = MachineMesh(msh.bandData.p_ag_virt, ...
                    msh.bandData.t_ag);
            end
        end
        
        function msh = generateMovingBand(msh, varargin)
            msh.bandData = AGtriangulation(msh, varargin{:});
        end
        
        function bl = isfield(msh, arg)
            %another compatibility fix for functions expecting a struct
            %with field rotel or statel.
            bl = false;
            if strcmp(arg, 'symmetrySectors')
                bl = true;
            elseif any(msh.namedElements.get(arg))
                bl = true;
            end
        end
        
        function msh = setMachineBoundaryNodes(msh, dims, varargin)
            %setMachineBoundaryNodes find and store boundary nodes.
            %
            % setMachineBoundaryNodes(msh, dims) is equivalent to calling
            % both
            % setMachineDirichletNodes(dims) and
            % setMachinePeriodicNodes()
            % in that order.
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            msh.setMachineDirichletNodes(dims, varargin{:});
            msh.setMachinePeriodicNodes(varargin{:});
        end
        
        function msh = setMachineDirichletNodes(msh, dims, varargin)
            %setMachineDirichletNodes Dirichlet nodes on outer and inner
            %boundaries.
            %
            % setMachineDirichletNodes(dims) and
            % setMachineDirichletNodes(dims, tol)
            % try to find the nodes on the outer (and possibly inner)
            % boundaries of the machines, based on the dimensions dims.D_so
            % (stator outer diameter) and optionally dims.D_ri (rotor inner diameter).
            % A defaults tolerance of 1e-4 is used,
            % unless a different one is given as input.
            %
            % Nodes are stored as namedNodes('Dirichlet').
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            
            TOL = 1e-4;
            if numel(varargin)
                TOL = varargin{1};
            end
            ro = dims.D_so / 2;
            try ri = dims.D_ri / 2; catch; ri = 0;  end
            n_out = find( abs(sum(msh.p.^2,1) - ro^2) < (TOL) );
            n_in = find( abs(sum(msh.p.^2,1) - ri^2) < (0.1*TOL) );
            
            n_dir = [sortSegmentEdges(msh.p, toRow(n_out)) ...
                sortSegmentEdges(msh.p, toRow(n_in))];
            
            msh.namedNodes.add('Dirichlet', n_dir);
        end
        
        function msh = setMachinePeriodicNodes(msh, varargin)
            %setMachinePeriodicNodes nodes on periodic boundaries.
            %
            % setMachinePeriodicNodes() and
            % setMachinePeriodicNodes(tol) try to find the nodes
            % belonging to the periodic boundaries, determined based on the
            % number of symmetry sectors. A default tolerance of 1e-4 is
            % used unless a different one is supplied as input.
            %
            % Master boundary is assumed to lie on the positive x-axis, and
            % its nodes are stored in namedNodes('Periodic_master').
            %
            % Slave boundary nodes are stored in
            % namedNodes('Periodic_slave')
            %
            % Dirichlet boundary nodes belong to neither, and must be set
            % before calling this function.
            %
            % Note that the boundaries are assumed SYMMETRIC by default,
            % i.e. that master node coordinates can be obtained from slave
            % node coordinates by a rotation around the origin. If this is
            % not the case (nodal coordinates don't match and/or the number
            % of nodes differs),
            % CALL msh.info.set('NonSymmetricBND', true) before running any
            % other function requiring boundary information.
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            
            TOL = 1e-4;
            if numel(varargin)
                TOL = varargin{1};
            end
            if msh.symmetrySectors == 1
                return;
            end
            n_master_cand = find( (msh.p(2,:) < TOL) & (msh.p(1,:)>0) );
            %n_master_cand = setdiff(n_master_cand, msh.namedNodes.get('Dirichlet'));
            
            %sorting based on radius and adding
            r2 = sum(msh.p(:,n_master_cand).^2,1); [~,I] = sort(r2);
            msh.namedNodes.add('Periodic_master', n_master_cand(I));
            
            if msh.symmetrySectors == 2
                n_slave_cand = find( (msh.p(1,:)<0) & (msh.p(2,:)<TOL) );
            else
                sectorAngle = 2*pi / msh.symmetrySectors;
                temp = cos(sectorAngle)*msh.p(2,:) - sin(sectorAngle)*msh.p(1,:);
                n_slave_cand = find(abs(temp)<TOL);
            end
            %n_slave_cand = setdiff(n_slave_cand, ...
            %    [msh.namedNodes.get('Dirichlet') msh.namedNodes.get('Periodic_master')]);
            r2 = sum(msh.p(:,n_slave_cand).^2,1);  [~,I] = sort(r2);
            
            msh.namedNodes.add('Periodic_slave', n_slave_cand(I));
        end
        
        function m2 = copy(msh, varargin)
            %copy A deep copy.
            %
            % See MeshBase.copy for documentation.
            %
            % (c) 2017 Antti Lehikoinen / Aalto University
            
            if numel(varargin)
                m2 = varargin{1};
            else
                m2 = MachineMesh();
            end
            m2 = copy@MeshBase(msh, m2);
            m2.symmetrySectors = msh.symmetrySectors;
            m2.periodicityCoeff = msh.periodicityCoeff;
            
            try
                names = fieldnames(msh.bandData);
                m2.bandData = struct();
                for k = 1:numel(names)
                    m2.bandData.(names{k}) = msh.bandData.(names{k});
                end
            catch
                %no bandData defined --> pass
            end
        end
        
        function this = to2ndOrder(this)
            %to2ndOrder Transform mesh to second order.
            %
            %The method transforms the mesh elements into non-curved second-order
            %elements. For the named nodes (airgap bnd, periodic, Dirichlet) to be
            %updated correctly, the following criteria must be met:
            %   - The method is only called after the nodes have been set with e.g.
            %       msh.namedNodes.set('n_ag_s', stator_ag_nodes);
            %   - The named nodes are ordered either radially (periodic boundaries) or
            %       circumferentially (airgap nodes, Dirichlet nodes excluding possible
            %       center node).
            %   - The lists of periodic nodes now contain all nodes on the periodic
            %       boundary. In other words, something like
            %           n_cl_s = setdiff(n_cl_s, n_dir_s);
            %       must NOT be called.
            %   - Air gap mesh generation (msh.generateMovingBand()) is only called
            %   AFTER msh.2ndOrder().
            
            Np = size(this.p, 2);
            [this.p, this.t] = to2ndOrder(this.p, this.t, this.edges, this.t2e);
            this.edges = [this.edges; (Np+1):(Np+size(this.edges,2))];
            
            this = msh_updateNamedNodes(this);
            
            this.elementType = Elements.triangle2;
        end
        
        function Sag = get_AGmatrix(this, rotorAngle, varargin)
            %air-gap matrix
            % WARNING: air-gap surfaces are assumed non-curved
            
            if isobject(this.bandData)
                Sag = this.bandData.get_AGmatrix(rotorAngle, varargin{:});
                return;
            end
            
            if isempty(this.bandData)
                if numel(varargin)
                    Sag = sparse(varargin{1}, varargin{1});
                else
                    Np = size(this.p, 2);
                    Sag = sparse(Np, Np);
                end
                return;
            end
            
            if this.elementType == Elements.triangle
                Sag = get_MovingBandMatrix(rotorAngle, this, varargin{:});
            elseif this.elementType == Elements.triangle2 || this.elementType == Elements.triangle2I
                %number of nodes to skip on rotor surface
                nodeShift = -floor( (rotorAngle - 1*this.bandData.shiftTol/2) / this.bandData.shiftTol ) - 1;
                nodeShift = 2*nodeShift;
                
                %shifting indices in the sorted list
                newPositions = mod(this.bandData.originalPositions_rotor - 1 + nodeShift, this.bandData.N_ag_r ) + 1;
                
                t_ag = this.bandData.t_ag;
                t_ag(this.bandData.inds_r) = this.bandData.sortedNodes_rotor(newPositions);
                
                p = this.bandData.p_ag_virt;
                p(:, this.bandData.sortedNodes_rotor) = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)] * ...
                    p(:, this.bandData.sortedNodes_rotor);
                
                %updating edge center node positions
                p(:, this.bandData.msh_ag.edges(3,:)) = ...
                    0.5*p(:, this.bandData.msh_ag.edges(1,:)) + ...
                    0.5*p(:, this.bandData.msh_ag.edges(2,:));
                this.bandData.msh_ag.p = p;
                this.bandData.msh_ag.t = t_ag;
                
                %figure(10); clf;
                %msh_triplot(this.bandData.msh_ag, [], 'b');
                
                Sag_c = MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), 1/(pi*4e-7), [], this.bandData.msh_ag);
                inds = 1:Sag_c.Nvals;
                Sag_c.E(inds) = Sag_c.E(inds) .* this.bandData.el_table(3, Sag_c.I(inds));
                Sag_c.I(inds) = this.bandData.el_table(2, Sag_c.I(inds));
                Sag_c.E(inds) = Sag_c.E(inds) .* this.bandData.el_table(3, Sag_c.J(inds));
                Sag_c.J(inds) = this.bandData.el_table(2, Sag_c.J(inds));
                if numel(varargin)
                    Sag = Sag_c.finalize(varargin{1}, varargin{1});
                else
                    Np = size(this.p, 2);
                    Sag = Sag_c.finalize(Np, Np);
                end
            end
        end
    end
end

