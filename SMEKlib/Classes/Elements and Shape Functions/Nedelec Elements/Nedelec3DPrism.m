classdef Nedelec3DPrism < handle
    %Nedelec3D class for 3D Nedelec shape functions.
    % 
    % Only supports tets, so far.
    % 
    % (c) 2018 Antti Lehikoinen / Aalto University
    properties
        op
    end
    
    methods
        function this = Nedelec3DPrism(varargin)
            if numel(varargin)
                opr = varargin{1};
            else
                opr = Operators.I;
            end
            this.op = opr;
        end
        
        function wref = eval_ref(this, k, x, msh)
            %eval_ref Evaluate in reference element.
            % 
            % Call syntax
            % wref = this.eval_ref(k, xref, ~), where
            %   k = number of edge, in [1,6]
            %   xref = coordinate in reference tetrahedron.
            %
            % Refence: "Higher Order Interpolatory Vector Bases on Prism
            % Elements"
            
            %barycentric coordinates etc
            p1 = x(2); p2 = 1 - x(1) - x(2); p3 = x(1); p4 = x(3); p5 = 1-x(3);
            d1 = [0; 1; 0]; d2 = [-1;-1; 0]; d3 = [1;0; 0]; d4 = [0;0;1]; d5 = [0;0;-1];
            
            l12 = [0; 0; 1]; l23 = l12; l13 = -l12;
            s12 = 1; s23 = 1; s13 = -1;
            
            l14 = [1; 0; 0]; l15 = -l14;
            s14 = 1; s15 = -1;
            
            l24 = [-1; 1; 0]; l25 = -l24;
            s24 = 1; s25 = -1;
            
            l34 = [0;-1;0]; l35 = -l34;
            s34 = 1; s35 = -1;
            
            J = 1;
            
            if this.op == Operators.I
                if k == 1
                    %14
                    wref = p5 * (p2*d3 - p3*d2);
                    wref = wref*s14;
                elseif k == 4
                    %15
                    wref = p4 * (p3*d2 - p2*d3);
                    wref = wref*s15;
                elseif k == 2
                    %24
                    wref = p5 * (p3*d1 - p1*d3);
                    wref = wref*s24;
                elseif k == 5
                    %25
                    wref = p4 * (p1*d3 - p3*d1);
                    wref = wref*s25;
                elseif k == 3
                    %34
                    wref = p5 * (p1*d2 - p2*d1);
                    wref = wref*s34;
                elseif k == 6
                    %35
                    wref = p4 * (p2*d1 - p1*d2);
                    wref = wref*s35;
                elseif k == 7
                    %13
                    wref = p2*d5;
                    wref = wref*s13;
                elseif k == 8
                    %12
                    wref = p3*d4;
                    wref = wref*s12;
                elseif k == 9
                    %23
                    wref = p1*d4;
                    wref = wref*s23;
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif this.op == Operators.curl
                if k == 1
                    %14
                    wref = p3*l25 - p2*l35 + 2*p5*l23;
                    wref = s14*wref;
                elseif k == 4
                    %15
                    wref = p2*l34 - p3*l24 - 2*p4*l23;
                    wref = s15*wref;
                elseif k == 2
                    %24
                    wref = p1*l35 - p3*l15 - 2*p5*l13;
                    wref = s24*wref;
                elseif k == 5
                    %25
                    wref = p3*l14 - p1*l34 + 2*p4*l13;
                    wref = s25*wref;
                elseif k == 3
                    %34
                    wref = p2*l15 - p1*l25 + 2*p5*l12;
                    wref = s34*wref;
                elseif k == 6
                    %35
                    wref = p1*l24 - p2*l14 - 2*p4*l12;
                    wref = s35*wref;
                elseif k == 7
                    %13
                    wref = l25;
                    wref = s13*wref;
                elseif k == 8
                    %12
                    wref = l34;
                    wref = s12*wref;
                elseif k == 9
                    %23
                    wref = l14;
                    wref = s23*wref;
                end
                wref = wref / J;
            end
        
        end
        
        function w = eval(this, k, X, msh, varargin)
            %eval Evaluate global shape function.
            % 
            % Call syntax 
            % w = this.eval(k, xref, msh, elements) or
            % w = this.eval(k, xref, msh, F, detF, elements)
            
            %getting mapping
            if size(varargin{1}, 1) == 9
                F = varargin{1};
                detF = varargin{2};
                elements = varargin{3};
            else
                elements = varargin{1};
                F = msh.getMappingMatrix(elements, X);
                detF = matrixDeterminant(F);
            end
            if isempty(elements)
                elements = 1:size(msh.t,2);
            end
            %elements        
            wref = this.eval_ref(k, X, msh);
            if this.op == Operators.I
                w = matrixTimesVector(F, wref, true, true, detF);
            elseif this.op == Operators.curl
                w = bsxfun(@times, matrixTimesVector(F, wref, false, false), 1./detF);
            end
            w = bsxfun(@times, w, sign(msh.elements2edges(k, elements)));
            
            %w = NedelecShapeFunction.evaluate(this.op, ShapeFunctions.nedelec, ...
            %    k, msh, X, elements, F, detF);
            
            %correcting orientation (already corrected inside
            %NedelecShapeFunction.evaluate)
            %w = bsxfun(@times, w, sign(msh.elements2edges(k, elements)));
        end
        
        function w = eval2(this, k, X, msh, varargin)
            %eval Evaluate global shape function.
            % 
            % Call syntax 
            % w = this.eval(k, xref, msh, elements) or
            % w = this.eval(k, xref, msh, F, detF, elements)
            
            %getting mapping
            if size(varargin{1}, 1) == 9
                F = varargin{1};
                detF = varargin{2};
                elements = varargin{3};
            else
                elements = varargin{1};
                F = msh.getMappingMatrix(elements, X);
                detF = matrixDeterminant(F);
            end
            if isempty(elements)
                elements = 1:size(msh.t,2);
            end
            
            w = NedelecShapeFunction.evaluate(this.op, ShapeFunctions.nedelec, ...
                k, msh, X, elements, F, detF);
            
            %correcting orientation (already corrected inside
            %NedelecShapeFunction.evaluate)
            %w = bsxfun(@times, w, sign(msh.elements2edges(k, elements)));
        end
        
        function [Nf, order, Nvars] = getData(this, msh)
            %getData get shape function and element data.
            %
            % Call syntax
            % [Nf, order, Nvars] = this.getData(msh), where
            %   Nf = number of shape functions per reference element
            %   order = order of shape functions
            %   Nvars = number of dofs in entire mesh.
            
            if this.op == Operators.I
                order = 1;
            else
                order = 0;
            end
            Nf = 9;
            Nvars = size(msh.edges, 2);
        end
        
        function inds = getIndices(~, k, msh, varargin)
            %getIndices Get indices of DoFs.
            %
            % Call syntax
            % inds = this.getIndices(k, msh)
            % inds = this.getIndices(k, msh, elements), where
            %   k = number of reference shape function.
            
            if ~numel(varargin) || ~any(varargin{1}) || (varargin{1}(1)<0)
                inds = abs(msh.elements2edges(k,:));
            else
                inds = abs(msh.elements2edges(k, varargin{1}));
            end
        end
    end
end