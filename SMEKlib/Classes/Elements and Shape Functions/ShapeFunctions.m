classdef ShapeFunctions < handle
    enumeration
        %shape function types
        nodal, nedelec,
        nedelec_all,
        nodal_all,
    end
    methods(Static)
        function [bt, Nf] = getBasisData(element, shapefunction)
           if element == Elements.triangle
               if shapefunction == ShapeFunctions.nodal
                   bt = BasisTypes.linear2; Nf = 3; return;
               else
                   error('Not yet implemented')
               end
           elseif element == Elements.prism
               if shapefunction == ShapeFunctions.nedelec
                   bt = BasisTypes.bilinearz3; Nf=9; return
               else
                   error('Not yet implemented')
               end
           else
               error('Element type not yet implemented');
           end
        end
        
        function v = evaluate(op, fun, kfun, msh, varargin)
            switch ShapeFunctions.getSuperClass(fun)
                case ShapeFunctions.nedelec_all
                    v = NedelecShapeFunction.evaluate(op, fun, kfun, msh, varargin{:});
            end
        end
        
        function C = getCoefficients(op, fun, msh)
            switch fun
                case ShapeFunctions.nedelec
                    switch msh.elementType
                        case Elements.prism
                            C = NedelecShapeFunction.PrismCoefficients(op);
                    end
            end
        end
        
        function c = getSuperClass(fun)
            if (fun == ShapeFunctions.nedelec)
                c = ShapeFunctions.nedelec_all; return;
            end
        end
                    
    end
end