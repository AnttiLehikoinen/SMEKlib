classdef Nodal2D < handle
    properties
        op
    end
    
    methods
        function this = Nodal2D(varargin)
            if numel(varargin)
                opr = varargin{1};
            else
                opr = Operators.I;
            end
            this.op = opr;
        end
		
        function Nref = eval_ref(this, k, x, msh)
            % evaluates the reference shape function (either identity or
            % gradient)
            if msh.elementType == Elements.triangle
                if this.op == Operators.I
                    coeffs = {[1 -1 -1]; [0 1 0]; [0 0 1]};
                    Nref = coeffs{k} * [ones(1, size(x,2)); x];
                    return
                else
                    grads = {[-1;-1]; [1;0]; [0; 1]};
                    Nref = grads{k};
                end
            elseif (msh.elementType == Elements.triangle2) || (msh.elementType == Elements.triangle2I)
                bx = [ones(1, size(x,2)); x;
                    x(1,:).^2; x(1,:).*x(2,:); x(2,:).^2];
                if this.op == Operators.I
                    coeffs = [     1    -3    -3     2     4     2
                        0    -1     0     2     0     0
                        0     0    -1     0     0     2
                        0     4     0    -4    -4     0
                        0     0     0     0     4     0
                        0     0     4     0    -4    -4];
                    Nref = coeffs(k,:) * bx;
                    return;
                else
                    cfx = [    -3     4     4     0     0     0
                        -1     4     0     0     0     0
                        0     0     0     0     0     0
                        4    -8    -4     0     0     0
                        0     0     4     0     0     0
                        0     0    -4     0     0     0];
                    cfy = [    -3     4     4     0     0     0
                        0     0     0     0     0     0
                        -1     0     4     0     0     0
                        0    -4     0     0     0     0
                        0     4     0     0     0     0
                        4    -4    -8     0     0     0];
                    Nref = [cfx(k,:); cfy(k,:)]*bx;
                    return
                end
            else
                error('Shape function not implemented.');
            end
        end
        
        function N = eval(this, k, X, msh, varargin)
            %eval evaluates the global shape function
            % call syntax 
            % eval(k, X, msh, F, detF) or
            % eval(k, X, msh, elements)
            
            Nref = this.eval_ref(k, X, msh);
            if this.op == Operators.I
                % identity needed? --> return
                N = Nref; return;
            end
            
            %getting mapping
            if size(varargin{1}, 1) == 4
                F = varargin{1};
                detF = varargin{2};
            else
                F = msh.getMappingMatrix(varargin{1}, X);
                detF = matrixDeterminant(F);
            end
            
            %evaluating global gradient / curl
            N = matrixTimesVector(F, Nref, true, true, detF); %gradient
            if this.op == Operators.grad
                return;
            elseif this.op == Operators.curl
                N = [0 1;-1 0] * N;
                return;
            else
                error('Invalid operator');
            end
        end
        
        function [Nf, order, Nvars] = getData(this, msh)
           %getData
           if msh.elementType == Elements.triangle
               Nf = 3;
               if this.op == Operators.I
                   order = 1;
               else
                   order = 0;
               end
               Nvars = size(msh.p, 2);
           elseif msh.elementType == Elements.triangle2 || msh.elementType == Elements.triangle2I
               Nf = 6;
               if this.op == Operators.I
                   order = 2;
               else
                   order = 1;
               end
               Nvars = size(msh.p, 2);
           else
                error(['Element type ' char(msh.elementType) ' not yet implemented.']);
            end
        end
        
        function inds = getIndices(~, k, msh, varargin)
            if ~numel(varargin) || ~any(varargin{1}) || (varargin{1}(1)<0)
                inds = msh.t(k,:);
            else
                inds = msh.t(k, varargin{1});
            end
        end
        
    end
end
    