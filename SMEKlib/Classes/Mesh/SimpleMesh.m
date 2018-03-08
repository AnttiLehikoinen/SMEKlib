classdef SimpleMesh < handle
    %SimpleMesh a mesh class reduced functionality.
    properties
        elementType
        p, t
    end
    
    methods
        function msh = SimpleMesh(p, t)
            msh.p = p;
            msh.t = t;
            switch size(t, 1)
                case 3
                    msh.elementType = Elements.triangle;
                case 6
                    msh.elementType = Elements.triangle2I;
                otherwise
                    error('Element type not included here yet.');
            end
        end
        
        function [F, F0] = getMappingMatrix(this, varargin)
            [F, F0] = msh_getMappingMatrix(this, varargin{:});
        end
        
    end
end