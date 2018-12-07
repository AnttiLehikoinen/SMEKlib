classdef ExtrudedPrismMesh < PrismMeshBase3D
    %ExtrudedMachineMesh A base class for a MachineMesh extruded into 3D.
    % 
    % (c) 2017 Antti Lehikoinen / Aalto University
    
    properties
        Nfaces_tri, Nfaces_square
        faces_tri, faces_square
    end
    
    methods
        function msh3 = ExtrudedPrismMesh(msh2, zs)
            %extrudedMachineMesh constructor.
            % 
            % Call syntax
            % msh3 = extrudedMachineMesh(msh2) to extrude a MachineMesh and
            % generate a new ExtrudedMachineMesh object msh3.
            
            Ne2 = size(msh2.t, 2);
            Np2 = size(msh2.p, 2);
            Nedges2 = size(msh2.edges, 2);
            Nl = numel(zs);
            
            Ne3 = Ne2 * (Nl-1);
            Np3 = Np2 * Nl;
            Nedges3 = Nedges2*Nl + (Nl-1)*Np2;
            Nfaces_square = Nedges2 * (Nl-1);
            Nfaces_tri = Ne2*Nl;
            
            msh3 = msh3@PrismMeshBase3D();
            msh3.elementType = Elements.prism;
            msh3.Nfaces_tri = Nfaces_tri;
            msh3.Nfaces_square = Nfaces_square;
            
            %setting nodes
            msh3.nodes = [repmat(msh2.p, 1, Nl);
                kron(zs, ones(1, Np2))];
            
            %setting elements
            msh3.elements = zeros(6, Ne3);
            for kl = 1:(Nl-1)
                msh3.elements(:, (1:Ne2) + (kl-1)*Ne2) = ...
                    [msh2.t + (kl-1)*Np2;
                    msh2.t + kl*Np2];
            end
            
            %setting edges
            msh3.edges = zeros(2, Nedges3);
            msh3.edges(:, 1:(Nl*Nedges2)) = kron(0:(Nl-1), Np2*ones(2, Nedges2)) + ...
                repmat(msh2.edges, 1, Nl);
            for kl = 1:(Nl-1)
                msh3.edges(:, (1:Np2) + Nl*Nedges2 + (kl-1)*Np2) = ...
                    [(1:Np2) + (kl-1)*Np2;
                    (1:Np2) + kl*Np2];
            end
            
            %setting faces
            t2e = order_t2e(msh2);
            msh3.faces_square = zeros(4, Nfaces_square);
            for kl = 1:(Nl-1)
                msh3.faces_square(:, (1:Nedges2) + (kl-1)*Nedges2) = ...
                    [1:Nedges2;
                    Nl*Nedges2 + (kl-1)*Np2 + msh2.edges(2,:);
                    -(kl*Nedges2 + (1:Nedges2));
                    -(Nl*Nedges2 + (kl-1)*Np2 + msh2.edges(1,:))];
            end
            msh3.faces_tri = zeros(3, Nfaces_tri);
            msh3.faces_tri(:, 1:Ne2) = t2e;
            for kl = 1:(Nl-1)
                msh3.faces_tri(:, (1:Ne2) + kl*Ne2) = ...
                    (kl*Nedges2 + abs(t2e)).*sign(t2e);
            end
            
            %setting faces-to-elements
            msh3.faces2elements = zeros(2, Nfaces_square + Nfaces_tri);
            msh3.faces2elements(:, 1:Nfaces_square) = (kron(0:(Nl-2), Nedges2*ones(2,Nedges2)) + ...
                repmat(abs(msh2.e2t), 1, Nl-1)) .* repmat(sign(msh2.e2t), 1, Nl-1);
            msh3.faces2elements(:, (Nfaces_square+1):end) = [kron(-1:(Nl-2), Ne2*ones(1,Ne2));
                kron(0:(Nl-1), Ne2*ones(1,Ne2))] + repmat(1:Ne2, 2, Nl);
            msh3.faces2elements(:, (Nfaces_square+1):(Nfaces_square+Ne2)) = ...
                flipud(msh3.faces2elements(:, (Nfaces_square+1):(Nfaces_square+Ne2)));
            msh3.faces2elements(2, (Nfaces_square+1):(Nfaces_square+Ne2)) = 0;
            msh3.faces2elements(2, (end-Ne2+1):end) = 0;
            
            %setting elements-to-edges
            msh3.elements2edges = zeros(9, Ne3);
            msh3.elements2edges(1:6, :) = ([kron(0:(Nl-2), Nedges2*ones(3, Ne2));
                kron(1:(Nl-1), Nedges2*ones(3, Ne2))] +  repmat(abs(t2e), 2, Nl-1)) .* ...
                repmat(sign(t2e), 2, Nl-1);
            msh3.elements2edges(7:9, :) = (Nl*Nedges2 + ...
                kron(0:(Nl-2), + Np2*ones(3, Ne2))) + ...
                repmat(msh2.t, 1, Nl-1);
            
            %setting elements-to-faces
            msh3.elements2faces = zeros(5, Ne3);
            for kl = 1:(Nl-1)
                msh3.elements2faces(1:3, (1:Ne2) + (kl-1)*Ne2) = ...
                    ((kl-1)*Nedges2+abs(t2e)).*sign(t2e);
            end
            msh3.elements2faces(4,:) = ...
                (kron(0:(Nl-2), Ne2*ones(1,Ne2)) + ...
                repmat(1:Ne2, 1, Nl-1)) + Nfaces_square;
            msh3.elements2faces(5,:) = ...
                (kron(1:(Nl-1), Ne2*ones(1,Ne2)) + ...
                repmat(1:Ne2, 1, Nl-1)) + Nfaces_square;
            
            % updating named quantities:
            %named faces from named edges
            keys = msh2.namedEdges.keys();
            for kk = 1:numel(keys)
                key = keys{kk};
                val = msh2.namedEdges.get(key);
                faces = kron(0:(Nl-2), Nedges2*ones(1, numel(val))) + ...
                    repmat(val, 1, Nl-1);
                msh3.namedFaces.add(key, faces);
            end
            %named elements
            keys = msh2.namedElements.keys();
            for kk = 1:numel(keys)
                try
                    key = keys{kk};
                    val = msh2.namedElements.get(key);
                    els = kron(0:(Nl-2), Ne2*ones(1, numel(val))) + ...
                        repmat(val, 1, Nl-1);
                    msh3.namedElements.add(key, els);
                catch
                    disp('Something failed with named elements.');
                end
            end
        end
        
        function y = t(msh3)
            y = msh3.elements;
        end
        
        function x0 = elementCenters(msh3, elem)
            if ~any(elem) || elem(1) < 0
                elem = 1:size(msh3.elements,2);
            end
            x0 = zeros(3, numel(elem));
            for k = 1:size(msh3.elements, 1)
                x0 = x0 + msh3.nodes(:, msh3.elements(k, elem));
            end
            x0 = x0 / size(msh3.elements, 1);
        end
        
        
        function [] = triplot(msh3, els, varargin)
            if els < 0
                els = 1:size(msh3.elements,2);
            end
            faces = msh3.elements2faces(1:3, els);
            edges = reshape(msh3.faces_square(:, abs(faces)), 1, []);
            edges = toRow(unique(abs(edges)));
            msh_plotEdges3D(msh3, edges, varargin{:});
        end
        
        function [F, varargout] = getMappingMatrix(this, elem, varargin)
            if ~any(elem) || elem(1) < 0
                elem = 1:size(this.elements,2);
            end
            F = zeros(9, numel(elem));
            F(1:2,:) = this.nodes(1:2, this.elements(2,elem)) - this.nodes(1:2, this.elements(1,elem));
            F(4:5,:) = this.nodes(1:2, this.elements(3,elem)) - this.nodes(1:2, this.elements(1,elem));
            F(9,:) = this.nodes(3, this.elements(4,elem)) - this.nodes(3, this.elements(1,elem));
            
            if nargout > 1
                varargout{1} = this.nodes(:, this.elements(1, elem));
            end
        end
    end
    
end