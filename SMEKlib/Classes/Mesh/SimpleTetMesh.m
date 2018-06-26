classdef SimpleTetMesh < handle
    %SimpleTetMesh minimal mesh of 3D tetrahedrons.
    % 
    % (c) 2018 Antti Lehikoinen / Aalto University
    
    properties
        p, t, e, t2e, elementType
    end
    
    methods
        function msh = SimpleTetMesh(p, t, e, t2e)
            %SimpleTetMesh mesh constructor.
            % 
            % Call syntax 
            % msh = SimpleTetMesh(p, t, e, t2e)
            msh.p = p;
            msh.t = t;
            msh.e = e;
            msh.t2e = t2e;
            msh.elementType = Elements.tet;
        end
        
        function [F, F0] = getMappingMatrix(this, elem, varargin)
            %getMappingMatrix Mapping matrix from reference to global
            %element.
            % 
            % Call syntax
            %   [F, F0] = getMappingMatrix(this, elem) OR
            %   [F, F0] = getMappingMatrix(this, elem, unused_args)
            
            if isempty(elem)
                elem = 1:size(this.t, 2);
            end
            
            F0 = this.p(:, this.t(1,elem));
            F = [this.p(:, this.t(2,elem)) - this.p(:, this.t(1,elem));
                this.p(:, this.t(3,elem)) - this.p(:, this.t(1,elem));
                this.p(:, this.t(4,elem)) - this.p(:, this.t(1,elem))];
        end
    end
end
                