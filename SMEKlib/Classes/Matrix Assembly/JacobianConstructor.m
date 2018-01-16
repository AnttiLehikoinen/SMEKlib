classdef JacobianConstructor < MatrixConstructorBase
    properties
        Fvals_test, Fvals_shape, DETF, Ires, msh, Ne
        symmetric, x_quad, w
    end
    methods
        function this = JacobianConstructor(msh, fun_test, fun_shape, symmetric)
            N_NONLINEARITY_ORDER = 2;

            this = this@MatrixConstructorBase(); %TODO: add size initialization
            
            %this.Np = size(msh.p, 2);
            this.Ne = size(msh.t, 2);
            this.msh = msh;
            
            %checking symmetricity
            this.symmetric = symmetric;
            
            %getting data for integration etc.
            [this.x_quad, this.w, Nrows, Ncols, N_test, N_shape] = ...
                get_AssemblyParameters(fun_test, fun_shape, msh, N_NONLINEARITY_ORDER);      
            this.Nrows = Nrows; this.Ncols = Ncols;
            N_quad = numel(this.w);
            
            %mapping and mapping determinant
            if Elements.isIsoparametric(msh.elementType)
                detF = zeros(N_quad, this.Ne);
                F = cell(N_quad, 1);
                for k_quad = 1:N_quad
                    F{k_quad} = msh.getMappingMatrix([], this.x_quad(:,k_quad));
                    detF(k_quad,:) = mappingDeterminant(F{k_quad});
                end
                detFun = @(kq)( detF(kq,:) );
                Ffun = @(kq)( F{kq} );
            else
                F = msh.getMappingMatrix();
                detF = matrixDeterminant(F);
                detFun = @(kq)( detF );
                Ffun = @(kq)( F );
            end
            this.DETF = abs(detF);
            
            %saving function values
            this.Fvals_test = cell(N_quad, N_test);
            for k_quad = 1:N_quad
                for k_test = 1:N_test
                    this.Fvals_test{k_quad, k_test} = fun_test.eval(k_test, this.x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad));
                end
            end
            if symmetric
                this.Fvals_shape = this.Fvals_test;
            else
                this.Fvals_shape = cell(N_quad, N_shape);
                for k_quad = 1:N_quad
                    for k_shape = 1:N_shape
                        this.Fvals_shape{k_quad,k_shape} = ...
                            fun_shape.eval(k_shape, this.x_quad(:,k_quad), msh, Ffun(k_quad), detFun(k_quad));
                    end
                end
            end
            
            %storing coordinates
            this.Ires = zeros(1, N_test*this.Ne);
            for k_shape = 1:N_shape
                if this.symmetric
                    %symmetric Jacobian
                    for k_test = k_shape:N_test
                        this.addCoordinates(fun_test.getIndices(k_test, msh), ...
                            fun_shape.getIndices(k_shape, msh));
                        if k_test ~= k_shape
                            this.addCoordinates(msh.t(k_test,:), msh.t(k_shape,:));
                        end
                    end
                else
                    %general Jacobian
                    for k_test = 1:N_test
                         this.addCoordinates(fun_shape.getIndices(k_shape, msh), ...
                             fun_test.getIndices(k_test, msh));
                    end
                end
            end
            for k_test = 1:N_test
                this.Ires( (k_test-1)*this.Ne + (1:this.Ne) ) = msh.t(k_test,:);
            end
            this.I = this.I(1:this.Ncoord);
            this.J = this.J(1:this.Ncoord);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [J, res] = eval(this, X, matfun)
            
            Ntest = size(this.Fvals_test, 2);
            Nshape = size(this.Fvals_shape, 2);
            
            %determinant function
            if size(this.DETF,1) == 1
                detFun = @(kq)( this.DETF );
            else
                detFun = @(kq)( this.DETF(kq,:) );
            end
            
            %computing entries
            E = zeros(1, Nshape*Ntest*this.Ne);
            Eres = zeros(1, Ntest*this.Ne);
            for k_quad = 1:numel(this.w)
                ri = 1;
                
                %evaluating B
                B = zeros(2, this.Ne);
                for k_shape = 1:Nshape
                    B = B + bsxfun(@times, this.Fvals_shape{k_quad, k_shape}, ...
                        transpose(X(this.msh.t(k_shape,:))) );
                end
                
                %evaluating material characteristics
                [nu, dnu] = matfun(B);
                if size(nu, 1) == 1
                    %computing H and dH/dB
                    H = bsxfun(@times, B, nu);

                    dH = [nu+2*dnu.*B(1,:).*B(1,:);
                        2*dnu.*B(1,:).*B(2,:);
                        2*dnu.*B(2,:).*B(1,:);
                        nu+2*dnu.*B(2,:).*B(2,:)];
                else
                    %function returns H and dH/dB
                    H = nu;
                    dH = dnu;
                end
                
                %accumulating matrix entries
                for k_shape = 1:Nshape
                    dhdb = matrixTimesVector(dH, this.Fvals_shape{k_quad, k_shape}, false, false);
                    if this.symmetric
                        for k_test = k_shape:Ntest
                            inds = (ri-1)*this.Ne + (1:this.Ne); ri = ri + 1;
                            Etemp = this.w(k_quad)*...
                                dotProduct(this.Fvals_test{k_quad, k_test}, dhdb) .* detFun(k_quad);
                            E(inds) = E(inds) + Etemp;
                            if k_test ~= k_shape
                                inds = (ri-1)*this.Ne + (1:this.Ne); ri = ri + 1;
                                E(inds) = E(inds) + Etemp;
                            end
                        end
                    else
                        for k_test = 1:Ntest
                            inds = (ri-1)*this.Ne + (1:this.Ne); ri = ri + 1;
                             E(inds) = E(inds) + this.w(k_quad)*...
                                dotProduct(this.Fvals_test{k_quad, k_test}, dhdb) .* detFun(k_quad);
                        end
                    end
                end
                
                %accumulating residual entries
                for k_test = 1:Ntest
                    inds_res = (k_test-1)*this.Ne + (1:this.Ne);
                    Eres(inds_res) = Eres(inds_res) + this.w(k_quad) * ...
                        dotProduct(this.Fvals_test{k_quad, k_test}, H) .* detFun(k_quad);
                end
            end
            
            %multiplication by determinant
            %E = E.*repmat(this.DETF, 1, Nshape*Ntest);
            %Eres = Eres.*repmat(this.DETF, 1, Ntest);
            
            %E(:, 1:10)
            
            %assembly
            N = size(X, 1);
            J = sparse(this.I, this.J, E, N, N);
            res = sparse(this.Ires, ones(1, numel(Eres)), Eres, N, 1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [Jrr, Jri, Jir, Jii, resr, resi] = eval_complex(this, X, matfun)
            Nshape = size(this.Fvals_shape, 2); Ntest = size(this.Fvals_test, 2);
            N = size(X,1)/2;
            
            %determinant function
            if size(this.DETF,1) == 1
                detFun = @(kq)( this.DETF );
            else
                detFun = @(kq)( this.DETF(kq,:) );
            end
            
            %computing entries
            E = zeros(4, Nshape*Ntest*this.Ne); Eres = zeros(2, Ntest*this.Ne);
            for k_quad = 1:numel(this.w)
                ri = 1;
                
                %evaluating B
                Br = zeros(2, this.Ne); Bi = zeros(2, this.Ne);
                for k_shape = 1:Nshape
                    Br = Br + bsxfun(@times, this.Fvals_shape{k_quad, k_shape}, ...
                        transpose(X(this.msh.t(k_shape,:))) );
                    Bi = Bi + bsxfun(@times, this.Fvals_shape{k_quad, k_shape}, ...
                        transpose(X(N + this.msh.t(k_shape,:))) );
                end
                
                %evaluating material characteristics
                [nu, dnu] = matfun( sum(Br.^2+Bi.^2,1) );
                
                Hr = bsxfun(@times, Br, nu); Hi = bsxfun(@times, Bi, nu);
                
                dRdR = [nu+2*dnu.*Br(1,:).*Br(1,:);2*dnu.*Br(1,:).*Br(2,:);
                    2*dnu.*Br(2,:).*Br(1,:); nu+2*dnu.*Br(2,:).*Br(2,:)];
                dRdI = [2*dnu.*Br(1,:).*Bi(1,:);2*dnu.*Br(1,:).*Bi(2,:);
                    2*dnu.*Br(2,:).*Bi(1,:); 2*dnu.*Br(2,:).*Bi(2,:)];
                dIdR = [2*dnu.*Bi(1,:).*Br(1,:);2*dnu.*Bi(1,:).*Br(2,:);
                    2*dnu.*Bi(2,:).*Br(1,:); 2*dnu.*Bi(2,:).*Br(2,:)];
                dIdI = [nu+2*dnu.*Bi(1,:).*Bi(1,:);2*dnu.*Bi(1,:).*Bi(2,:);
                    2*dnu.*Bi(2,:).*Bi(1,:); nu+2*dnu.*Bi(2,:).*Bi(2,:)];

                
                %accumulating matrix entries
                for k_shape = 1:Nshape
                    dRR = matrixTimesVector(dRdR, this.Fvals_shape{k_quad, k_shape}, false, false);
                    dRI = matrixTimesVector(dRdI, this.Fvals_shape{k_quad, k_shape}, false, false);
                    dIR = matrixTimesVector(dIdR, this.Fvals_shape{k_quad, k_shape}, false, false);
                    dII = matrixTimesVector(dIdI, this.Fvals_shape{k_quad, k_shape}, false, false);
                    for k_test = 1:Ntest
                        inds = (ri-1)*this.Ne + (1:this.Ne); ri = ri + 1;
                        E(1,inds) = E(1,inds) + this.w(k_quad)*...
                            dotProduct(this.Fvals_test{k_quad, k_test}, dRR) .* detFun(k_quad);
                        E(2,inds) = E(2,inds) + this.w(k_quad)*...
                            dotProduct(this.Fvals_test{k_quad, k_test}, dRI) .* detFun(k_quad);
                        E(3,inds) = E(3,inds) + this.w(k_quad)*...
                            dotProduct(this.Fvals_test{k_quad, k_test}, dIR) .* detFun(k_quad);
                        E(4,inds) = E(4,inds) + this.w(k_quad)*...
                            dotProduct(this.Fvals_test{k_quad, k_test}, dII) .* detFun(k_quad);
                    end
                end
                
                %accumulating residual entries
                for k_test = 1:Ntest
                    inds_res = (k_test-1)*this.Ne + (1:this.Ne);
                    Eres(1, inds_res) = Eres(1, inds_res) + this.w(k_quad) * ...
                        dotProduct(this.Fvals_test{k_quad, k_test}, Hr) .* detFun(k_quad);
                    Eres(2, inds_res) = Eres(2, inds_res) + this.w(k_quad) * ...
                        dotProduct(this.Fvals_test{k_quad, k_test}, Hi) .* detFun(k_quad);
                end
            end
            
            %multiplication by determinant
            %E = E.*repmat(this.DETF, 4, Nshape*Ntest);
            %Eres = Eres.*repmat(this.DETF, 2, Ntest);
            
            %assembly
            Jrr = sparse(this.I, this.J, E(1,:), N, N);
            Jri = sparse(this.I, this.J, E(2,:), N, N);
            Jir = sparse(this.I, this.J, E(3,:), N, N);
            Jii = sparse(this.I, this.J, E(4,:), N, N);
            resr = sparse(this.Ires, ones(1, size(Eres,2)), Eres(1,:), N, 1);
            resi = sparse(this.Ires, ones(1, size(Eres,2)), Eres(2,:), N, 1);
        end
    end
end