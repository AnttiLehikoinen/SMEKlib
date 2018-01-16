classdef NedelecShapeFunction < handle    
    methods (Static)
        function v = evaluate(op, fun, kfun, msh, X, els, varargin)
            if (fun == ShapeFunctions.nedelec) && (msh.elementType == Elements.prism)
                C = NedelecShapeFunction.PrismCoefficients(op);
                vloc = C{kfun}*...
                    PolynomialBasis3D.bteval(X, BasisTypes.bilinearz3);
                if ~any(els) || els(1) < -1
                    els = 1:size(msh.elements,2);
                end                
                if numel(varargin) == 2
                    F = varargin{1};
                    detF = varargin{2};
                else
                    F = msh.getMappingMatrix(els);
                    detF = matrixDeterminant(F);
                end
                if op == Operators.I
                    v = matrixTimesVector(F, vloc, true, true, detF);
                elseif op == Operators.curl
                    v = bsxfun(@rdivide, matrixTimesVector(F, vloc, false, false), detF);
                else
                    error(['Operator ' char(op) ' not yet implemented.']);
                end
                v = bsxfun(@times, v, sign(msh.elements2edges(kfun, els)));
            end
        end
                    
                    
        
        function C = PrismCoefficients(op)
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