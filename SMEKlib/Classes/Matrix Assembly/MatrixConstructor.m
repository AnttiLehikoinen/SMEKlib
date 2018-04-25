classdef MatrixConstructor < MatrixConstructorBase
    properties
        %for storing values
        funvals_test, funvals_shape, determinant
    end
    methods
        function this = MatrixConstructor(varargin)
            this = this@MatrixConstructorBase(); %TODO: add size initialization
            if numel(varargin)
                this = assemble_matrix(this, varargin{:});
            end
        end
        
        function this = assemble_vector(this, fun_test, k, v, elements, msh, varargin)
            this.assemble_matrix(fun_test, IDfun(k), v, elements, msh, varargin{:});
        end
        
        function this = assemble_matrix(this, fun_test, fun_shape, v, elements, msh, varargin)
            %assemble_matrix.
            % 
            % assemble_matrix(fun_test, fun_shape, v, elements, msh)
            
            %no elements listed --> going over all
            if ~any(elements)
                elements = 1:size(msh.t,2);
            end
            Ne = numel(elements);

            %size of v(x) too small --> assuming constant value
            if numel(v) == 1
                v = v(1) * ones(1, numel(elements));
            elseif size(v,2) ~= numel(elements)
                v = v(elements);
            end
            
            %checking symmetricity
            if strcmp(class(fun_test), class(fun_shape)) && (fun_test.op == fun_shape.op)
                symmetric = true;
            else
                symmetric = false;
            end
            
            %getting data for integration etc.
            [x_quad, w_quad, Nrows, Ncols, N_test, N_shape] = ...
                get_AssemblyParameters(fun_test, fun_shape, msh);
            N_quad = numel(w_quad);
            this.Nrows = Nrows; this.Ncols = Ncols;
            
            %mapping and mapping determinant
            if Elements.isIsoparametric(msh.elementType)
                detF = zeros(N_quad, Ne);
                F = cell(N_quad, 1);
                for k_quad = 1:N_quad
                    F{k_quad} = msh.getMappingMatrix(elements, x_quad(:,k_quad));
                    detF(k_quad,:) = mappingDeterminant(F{k_quad});
                end
                detFun = @(kq)( detF(kq,:) );
                Ffun = @(kq)( F{kq} );
            else
                F = msh.getMappingMatrix(elements);
                detF = matrixDeterminant(F);
                detFun = @(kq)( detF );
                Ffun = @(kq)( F );
            end
            
            %pre-computing shape functions
            Fvals_test = cell(N_quad, N_test);
            for k_quad = 1:N_quad
                for k_test = 1:N_test
                    Fvals_test{k_quad, k_test} = fun_test.eval(k_test, x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad), elements);
                end
            end
            if symmetric
                Fvals_shape = {};
            else
                Fvals_shape = cell(N_quad, N_shape);
                for k_quad = 1:N_quad
                    for k_shape = 1:N_shape
                        Fvals_shape{k_quad,k_shape} = ...
                            fun_shape.eval(k_shape, x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad), elements);
                    end
                end
            end
            
            %integrating
            rbias = this.Nvals;
            for k_quad = 1:N_quad
                ri = 0;
                for k_test = 1:N_test
                    if symmetric
                        for k_shape = k_test:N_shape                            
                            %Etemp = w_quad(k_quad)*dotProduct(Fvals_test{k_quad, k_test}, Fvals_test{k_quad, k_shape}).* v .* abs(detFun(k_quad));
                            Etemp = w_quad(k_quad)*FEdotProduct(Fvals_test{k_quad, k_test}, v, Fvals_test{k_quad, k_shape}) .* abs(detFun(k_quad));
                            this.addValues(Etemp, ri*Ne + rbias); ri = ri + 1;
                            if k_shape ~= k_test
                                this.addValues(Etemp, ri*Ne + rbias); ri = ri + 1;
                            end
                        end
                    else
                        for k_shape = 1:N_shape
                            %Etemp = w_quad(k_quad)*dotProduct(Fvals_test{k_quad, k_test}, Fvals_shape{k_quad, k_shape}).* v .* abs(detFun(k_quad));
                            Etemp = w_quad(k_quad)*FEdotProduct(Fvals_test{k_quad, k_test}, v, Fvals_shape{k_quad, k_shape}) .* abs(detFun(k_quad));
                            this.addValues(Etemp, ri*Ne + rbias); ri = ri + 1;
                            %error('Should not go here.')
                        end
                    end
                end
            end
            
            %adding indices
            for k_test = 1:N_test
                if symmetric
                    for k_shape = k_test:N_shape
                        this.addCoordinates( fun_test.getIndices(k_test, msh, elements), fun_shape.getIndices(k_shape, msh, elements) );
                        if k_shape ~= k_test
                            this.addCoordinates( fun_shape.getIndices(k_shape, msh, elements), fun_test.getIndices(k_test, msh, elements) );
                        end
                    end
                else
                    for k_shape = 1:N_shape
                        this.addCoordinates( fun_test.getIndices(k_test, msh, elements), fun_shape.getIndices(k_shape, msh, elements) );
                    end
                end
            end
            
        end
        
        function this = assemble_with_fun(this, fun_test, fun_shape, fun_v, elements, msh, order_v)
            %no elements listed --> going over all
            if ~any(elements)
                elements = 1:size(msh.t,2);
            end
            Ne = numel(elements);
            
           
            %getting data for integration etc.
            [x_quad, w_quad, Nrows, Ncols, N_test, N_shape] = ...
                get_AssemblyParameters(fun_test, fun_shape, msh, 0, order_v);
            N_quad = numel(w_quad);
            this.Nrows = Nrows; this.Ncols = Ncols;
            
            %mapping and mapping determinant
            if Elements.isIsoparametric(msh.elementType)
                detF = zeros(N_quad, Ne);
                F = cell(N_quad, 1);
                X_quad = cell(N_quad, 1);
                for k_quad = 1:N_quad
                    [F{k_quad}, X_quad{k_quad}] = msh.getMappingMatrix(elements, x_quad(:,k_quad));
                    detF(k_quad,:) = mappingDeterminant(F{k_quad});
                end
                detFun = @(kq)( detF(kq,:) );
                Ffun = @(kq)( F{kq} );
                Xfun = @(kq)( X_quad{kq} );
            else
                [F, F0] = msh.getMappingMatrix(elements);
                Xfun = @(kq)( F0 + matrixTimesVector(F, x_quad(:,kq), false, false) );
                detF = matrixDeterminant(F);
                detFun = @(kq)( detF );
                Ffun = @(kq)( F );
            end
            
            %computing test and shape function values
            Fvals_test = cell(N_quad, N_test);
            for k_quad = 1:N_quad
                for k_test = 1:N_test
                    Fvals_test{k_quad, k_test} = fun_test.eval(k_test, x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad));
                end
            end
            Fvals_shape = cell(N_quad, N_shape);
            for k_quad = 1:N_quad
                for k_shape = 1:N_shape
                    Fvals_shape{k_quad,k_shape} = ...
                        fun_shape.eval(k_shape, x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad));
                end
            end
            
            %integrating
            rbias = this.Nvals;
            for k_quad = 1:N_quad
                ri = 0;
                v = fun_v(x_quad(:,k_quad), Xfun(k_quad), Fvals_shape{k_quad, :});
                for k_test = 1:N_test
                    for k_shape = 1:N_shape
                        Etemp = w_quad(k_quad)*dotProduct(Fvals_test{k_quad, k_test}, Fvals_shape{k_quad, k_shape}).* v .* abs(detFun(k_quad));
                            this.addValues(Etemp, ri*Ne + rbias); ri = ri + 1;
                    end
                end                
            end
            
            %adding indices
            for k_test = 1:N_test
                for k_shape = 1:N_shape
                        this.addCoordinates( fun_test.getIndices(k_test, msh, elements), fun_shape.getIndices(k_shape, msh, elements) );
                end
            end
        end
    end
end