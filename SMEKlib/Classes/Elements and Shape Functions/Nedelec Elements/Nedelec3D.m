classdef Nedelec3D < handle
    %Nedelec3D class for 3D Nedelec shape functions.
    % 
    % Only supports tets, so far.
    % 
    % (c) 2018 Antti Lehikoinen / Aalto University
    properties
        op
    end
    
    methods
        function this = Nedelec3D(varargin)
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
            % Refence: "Fast MATLAB assembly of FEM matrices in 2D and 3D:
            % Edge elements"
            
            if this.op == Operators.curl
                switch k
                    case 1
                        wref = [0; -2; 1+1];
                    case 2
                        wref = [1+1; 0; -1-1];
                    case 3
                        wref = [-1-1; --1+1; 0];
                    case 4
                        wref = [0; 0; 1--1];
                    case 5
                        wref = [1--1; 0; 0];
                    case 6
                        wref = [0; --1+1; 0];
                end
                return
            end           
            
            x1 = x(1,:); x2 = x(2,:); x3 = x(3,:);
            switch k
                case 1
                    wref = [1-x3-x2; x1; x1];
                case 2
                    wref = [x2; 1-x3-x1; x2];
                case 3
                    wref = [x3; x3; 1-x2-x1];
                case 4
                    wref = [-x2; x1; 0*x1];
                case 5
                    wref = [0*x1; -x3; x2];
                case 6
                    wref = [x3; 0*x1; -x1];
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
            
            %correcting orientation
            w = bsxfun(@times, w, sign(msh.t2e(k, elements)));
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
            Nf = 6;
            Nvars = size(msh.e, 2);
        end
        
        function inds = getIndices(~, k, msh, varargin)
            %getIndices Get indices of DoFs.
            %
            % Call syntax
            % inds = this.getIndices(k, msh)
            % inds = this.getIndices(k, msh, elements), where
            %   k = number of reference shape function.
            
            if ~numel(varargin) || ~any(varargin{1}) || (varargin{1}(1)<0)
                inds = abs(msh.t2e(k,:));
            else
                inds = abs(msh.t2e(k, varargin{1}));
            end
        end
    end
end