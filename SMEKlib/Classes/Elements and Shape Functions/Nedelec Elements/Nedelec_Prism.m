classdef Nedelec_Prism < handle
    methods (Static)
        function C = getCoefficients(op)
            persistent CI Ccurl
            switch op
                case Operators.I
                    if isempty(CI)
                    	CI = preCompute_referenceShapeFunctions(Elements.prism, ShapeFunctions.nedelec);
                    end
                    C = CI; return;
                case Operators.curl
                    if isempty(Ccurl)
                        Ctemp = Nedelec_Prism.getCoefficients(Operators.I);
                        Ccurl = PolynomialBasis3D.curl_coeffs(Ctemp, Elements.prism, ShapeFunctions.nedelec);
                    end
                    C = Ccurl; return
                otherwise
                    error('Operator type not yet implemented')
            end

        end
    end
end