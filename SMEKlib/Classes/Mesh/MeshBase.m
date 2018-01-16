classdef MeshBase < handle
    %MeshBase base class for mesh objects.
    %
    % Has the following properties:
    %   p = nodal coordinates as a 2xNp array
    %   t = elements as a 3xNe array
    %   edges = ordered edges as a 2xNedges array, with edges(1,k)<edges(2,k)
    %   t2e = 3xNe incidence array: listing the 3 edges belonging to each element
    %   e2t = 2xNedges incidence array: from edges to elements. Second is zero for boundary edges
    %   matel = element material indices (defaults to zero)
    %
    %   namedNodes = SLContainer of named nodes (optional)
    %   namedElements = SLContainer of named elements (optional)
    %   namedEdges = SLContainer of named edges (optional)
    %
    %   info = a general SLContainer for key-value pairs. 
    %
    % (c) 2017 Antti Lehikoinen / Aalto University
    properties
        elementType
        p, t
        edges, t2e, e2t
        matel
        
        namedNodes, namedElements, namedEdges, info
    end
    
    methods
        function msh = MeshBase(p, t, varargin)
            msh.namedNodes = SLContainer();
            msh.namedElements = SLContainer();
            msh.namedEdges = SLContainer();
            msh.info = SLContainer();
            if nargin == 0
                %no arguments given; pass
            elseif nargin >= 2
                %initialization function; see help inittri for equivalent
                %documentation
                if size(t,1)~=3
                    %t = t';
                end
                %t = sort(t, 1);

                %defining edges an' pals
                [medges, me2t, mt2e] = getEdges(t);
                msh.p = p; msh.t = t; msh.edges = medges; msh.e2t = me2t;
                msh.t2e = mt2e;
                
                msh.matel = zeros(1, size(msh.t, 2)); %only air
            end
            if nargin == 3
                msh.matel = varargin{1};
            end
        end
        
        function msh = setDirichletNodes(msh)
            %automatic determination of Dirichlet nodes
            n_dir = unique( msh.edges(:, ~msh.e2t(2,:)) );
            msh = addNamedNodes(msh, 'Dirichlet', n_dir);
        end
        
        function msh2 = copy(msh, varargin)
            %copy A deep copy.
            %
            % msh.copy() returns a deep copy.
            % msh.copy(msh2) copies the data from msh to msh2. Mainly used
            % in subclasses of MeshBase.
            %
            % (c) 2017 Antti Lehikoinen

            if numel(varargin)
                msh2 = varargin{1};
            else
                msh2 = MeshBase();
            end
            msh2.p = msh.p;
            msh2.t = msh.t;
            msh2.edges = msh.edges;
            msh2.e2t = msh.e2t;
            msh2.t2e = msh.t2e;
            
            msh2.namedNodes = msh.namedNodes.copy();
            msh2.namedElements = msh.namedElements.copy();
            msh2.namedEdges = msh.namedEdges.copy();
            msh2.info = msh.info.copy();
        end
        
        function msh = split_edges(msh, edges)
            msh_split_edges2(msh, edges);
        end
        
        function [F, F0] = getMappingMatrix(this, varargin)
            if this.elementType == Elements.triangle || this.elementType == Elements.triangle2
                %getMappingMatrix A mesh-specific mapping matrix.            
                [F, F0] = mappingTerms(this, varargin{:});
            elseif this.elementType == Elements.triangle2I
                % [dF, F0] = getMappingMatrix(elements, x0)
                %isoparametric mapping F; returning dF/dxref(x0) and F(x0)
                Nref = Nodal2D(); gradN = Nodal2D(Operators.grad);
                if isempty(varargin{1})
                    elements = 1:size(this.t,2);
                else
                    elements = varargin{1};
                end
                x0 = varargin{2};
                F = zeros(4, numel(elements)); F0 = zeros(2, numel(elements));
                for kf = 1:size(this.t,1)
                    ngrad = gradN.eval_ref(kf, x0, this);
                    F = F + [this.p(:, this.t(kf,elements))*ngrad(1);
                        this.p(:, this.t(kf,elements))*ngrad(2)];
                    F0 = F0 + this.p(:, this.t(kf,elements))*Nref.eval_ref(kf, x0, this);
                end
            end
        end
    end
    
end