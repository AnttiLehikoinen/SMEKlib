classdef PolynomialBasis < handle
    properties
        xorder, yorder, zorder, basisType
    end
    %%{
    methods(Static)
        function bt = getBasisType(element, shapefunction)
           if element == Elements.triangle
               if shapefunction == ShapeFunctions.nodal
                   bt = PolynomialBasis.linear2; return;
               else
                   error('Not yet implemented')
               end
           elseif element == Elements.prism
               if shapefunction == ShapeFunctions.nedelec
                   bt = PolynomialBasis.bilinearz3; return
               else
                   error('Not yet implemented')
               end
           else
               error('Element type not yet implemented');
           end
        end
        function [xorder, yorder, zorder] = getOrdersFromBT(bt)
            switch bt
                case PolynomialBasis.bilinearz3
                    xorder = [0 1 0 0 1 0]';
                    yorder = [0 0 1 0 0 1]';
                    zorder = [0 0 0 1 1 1]';
                otherwise
                    error('Todo.');
            end
        end
        function [xorder, yorder, zorder] = getOrders(element, shapefunction)
            [xorder, yorder, zorder] = PolynomialBasis.getOrdersFromBT( ...
                PolynomialBasis.getBasisType(element, shapefunction) );
        end
           
        
        function v = eseval(X, element, shapefunction)
            [xorder, yorder, zorder] = PolynomialBasis.getOrders(element, shapefunction);            
            v = bsxfun(@power, X(1,:), xorder) .* ...
                bsxfun(@power, X(2,:), yorder) .* ...
                bsxfun(@power, X(3,:), zorder);
        end
        function v = bteval(X, bt)
            [xorder, yorder, zorder] = PolynomialBasis.getOrdersFromBT(bt);
            v = bsxfun(@power, X(1,:), xorder) .* ...
                bsxfun(@power, X(2,:), yorder) .* ...
                bsxfun(@power, X(3,:), zorder);
        end

        function C2 = curl_coeffs(C, varargin)
            if numel(varargin) == 1
                bt = varargin{1};
            else
                bt = PolynomialBasis.getBasisType(varargin{1}, varargin{2});
            end
            [xorder, yorder, zorder] = PolynomialBasis.getOrdersFromBT(bt);
            
            [dxM, dyM, dzM] = PolynomialBasis.derivated_coeffs(xorder, yorder, zorder);
            
            if isa(C, 'cell')
                C2 = cell(size(C));
                for k = 1:numel(C2)
                    C2{k} = [C{k}(3,:)*dyM - C{k}(2,:)*dzM;
                        -C{k}(3,:)*dxM + C{k}(1,:)*dzM;
                        C{k}(2,:)*dxM - C{k}(1,:)*dyM];
                end
            else
                C2 = [C(3,:)*dyM - C(2,:)*dzM;
                    -C(3,:)*dxM + C(1,:)*dzM;
                    C(2,:)*dxM - C(1,:)*dyM];
            end
        end
        
        function [dxM, dyM, dzM] = derivated_coeffs(xorder, yorder, zorder)
            cold = [xorder yorder zorder];
            M = eye(numel(xorder));
            
            dxord = max(xorder-1, 0);
            cnew = [dxord yorder zorder];
            [~, Locb] = ismember(cnew, cold, 'rows');
            dxM = M(Locb,:);
            dxM = bsxfun(@times, dxM, xorder);
            
            dord = max(yorder-1, 0);
            cnew = [xorder dord zorder];
            [~, Locb] = ismember(cnew, cold, 'rows');
            dyM = M(Locb,:);
            dyM = bsxfun(@times, dyM, yorder);
            
            dord = max(zorder-1, 0);
            cnew = [xorder yorder dord];
            [~, Locb] = ismember(cnew, cold, 'rows');
            dzM = M(Locb,:);
            dzM = bsxfun(@times, dzM, zorder);
        end
    end
    
    %}
    %non-static methods
    methods
        %%{
        function Pb = PolynomialBasis(varargin)
            
            
            if numel(varargin) == 1
                bt = varargin{1};
            else
                bt = PolynomialBasis.getBasisType(varargin{1}, varargin{2});
            end
            Pb.basisType = bt;
            [xorder, yorder, zorder] = PolynomialBasis.getOrdersFromBT(bt);
            Pb.xorder = xorder; Pb.yorder = yorder; Pb.zorder = zorder;
        end
        %}
        function v = eval(this, X)
            v = bsxfun(@power, X(1,:), this.xorder) .* ...
                bsxfun(@power, X(2,:), this.yorder) .* ...
                bsxfun(@power, X(3,:), this.zorder);
        end
    end

    enumeration
        linear2 (PolynomialBasis.linear2)
        bilinearz3 (PolynomialBasis.bilinearz3)      
    end
end             
        