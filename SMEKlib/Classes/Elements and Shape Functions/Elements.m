classdef Elements < handle
    enumeration
        %element types
        triangle, prism, composite,
        triangle2 %non-curved second-order
        triangle2I %isoparametric second-order
    end
    
    methods (Static)
        
        function bl = isIsoparametric(type)
            if type == Elements.triangle2I
                bl = true; return;
            end
            bl = false;
        end
        function bl = isTriangle(type)
            if type == Elements.triangle || type == Elements.triangle2 || type == Elements.triangle2I
                bl = true;
                return;
            end
            bl = false;
        end
        
        function X = refPoints_nodes(type)
            switch type
                case Elements.prism
                    X = [0 0 0; 1 0 0; 0 1 0;
                        0 0 1; 1 0 1; 0 1 1]';
                otherwise
                    error('');
            end
        end
        
        function X = refPoints_edges(type, varargin)
            switch type
                case Elements.prism
                    tau = [1 0 0; -1 1 0; 0 -1 0; 
                        1 0 0; -1 1 0; 0 -1 0;
                        0 0 1; 0 0 1; 0 0 1]';
                    tau0 = Elements.refPoints_nodes(Elements.prism);
                    Np = 7;
                    X = repmat(tau0(:,[1:end 1:3]), 1, Np) + kron(linspace(0,1, Np), tau);
                otherwise
                    error('');
            end
        end
    end
end