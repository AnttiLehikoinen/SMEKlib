function M_int = assemble_boundaryIntegralMatrix(op1, f1, n_bnd, elements, msh)


if ~(strcmp(op1, '') || strcmp(op1, 'grad'))
    error('Operator not implemented!');
elseif ~strcmp(f1, 'nodal')
    error('Shape function not implemented!');
end

Ne = size(n_bnd, 2);
Np = size(msh.p, 2);

%getting edge start- and end-points
X_start = msh.p(:, n_bnd(1,:));
X_end = msh.p(:, n_bnd(1,:));
%l_e = sum( (X_end-X_start).^2, 1).^0.5; %edge lengths (taken into account
%in the normal vector)

%determining mappings
[F, F0] = mappingTerms(msh_master, elements); 
detF = mappingDeterminant(F);

%getting segment start and end points in respective reference elements
xref_start = mappingTimesVector(X_start, 1, 0, F, -F0, detF);
xref_end = mappingTimesVector(X_start, 1, 0, F, -F0, detF);

%getting right-hand normal vectors
ne = [0 1;-1 0] * (X_start - X_end);

%defining function handles and indices
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] };
fh_test = @(ind_n, X)( phi_ref{ind_n}*X );
if strcmp(op1, 'grad')
    fh_shape = @(ind_n, X)( dotProduct(mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F, [], detF), ne) );
    NO_shape = 3;
    INDM_shape = msh.t(:,elements);
elseif strcmp(op1, '')
    fh_shape = fh_test;
    NO_shape = 2;
    INDM_shape = n_bnd;
else
    error('This should never happen...');
end
NO_test = 2;
INDM_test = n_bnd;

%getting reference 1D quad points
N_quad = 5;
[x_quad, W_quad] = internal_getQuadPoints(N_quad);
W_quad = W_quad/2; %scaling for interval length

%setting arrays
E = zeros(1, Ne*NO_test*NO_shape); I = E; J = E;
for k_quad = 1:N_quad
    x_quad = [(xref_end - xref_start)*x_quad(k_quad) / 2 + ...
        (xref_end + xref_start) / 2;
        ones(1, Ne)];
    
    ri = 1;
    for k_test = 1:NO_f
        for k_shape = 1:NO_f
            inds = ((ri-1)*Ne + 1):(ri*Ne);
            ri = ri + 1;
            
            E(inds) = E(inds) + W_quad(k_quad) * fh_test(k_test, x_quad) .* fh_shape(k_shape, x_quad);
        end
    end
end

%setting indices
ri = 1;
for k_test = 1:NO_test
    for k_shape = 1:NO_shape
        inds = ((ri-1)*Ne + 1):(ri*Ne);
        
        I(inds) = INDM_test(k_test,:);
        J(inds) = INDM_shape(k_shape,:);
        
        ri = ri + 1;
    end
end

M_int = sparse(I,J,E, Np, Np);

end

function [x_quad, W_quad] = internal_getQuadPoints(N_quad)

%returns 1D Gaussian quadrature points;
%copy-pasted from https://pomax.github.io/bezierinfo/legendre-gauss.html

switch N_quad
    case 2
        d = ...
            [1	1.0000000000000000	-0.5773502691896257
            2	1.0000000000000000	0.5773502691896257];
    case 3
        d = ...
            [1	0.8888888888888888	0.0000000000000000
            2	0.5555555555555556	-0.7745966692414834
            3	0.5555555555555556	0.7745966692414834];
    case 4
        d = ...
            [1	0.6521451548625461	-0.3399810435848563
            2	0.6521451548625461	0.3399810435848563
            3	0.3478548451374538	-0.8611363115940526
            4	0.3478548451374538	0.8611363115940526];
    case 5
        d = ...
            [1	0.5688888888888889	0.0000000000000000
            2	0.4786286704993665	-0.5384693101056831
            3	0.4786286704993665	0.5384693101056831
            4	0.2369268850561891	-0.9061798459386640
            5	0.2369268850561891	0.9061798459386640];
    case 6
        d = ...
            [1	0.3607615730481386	0.6612093864662645
            2	0.3607615730481386	-0.6612093864662645
            3	0.4679139345726910	-0.2386191860831969
            4	0.4679139345726910	0.2386191860831969
            5	0.1713244923791704	-0.9324695142031521
            6	0.1713244923791704	0.9324695142031521];
    case 7
        d = ...
            [1	0.4179591836734694	0.0000000000000000
            2	0.3818300505051189	0.4058451513773972
            3	0.3818300505051189	-0.4058451513773972
            4	0.2797053914892766	-0.7415311855993945
            5	0.2797053914892766	0.7415311855993945
            6	0.1294849661688697	-0.9491079123427585
            7	0.1294849661688697	0.9491079123427585];
    case 8
        d = ...
            [1	0.3626837833783620	-0.1834346424956498
            2	0.3626837833783620	0.1834346424956498
            3	0.3137066458778873	-0.5255324099163290
            4	0.3137066458778873	0.5255324099163290
            5	0.2223810344533745	-0.7966664774136267
            6	0.2223810344533745	0.7966664774136267
            7	0.1012285362903763	-0.9602898564975363
            8	0.1012285362903763	0.9602898564975363];
    case 9
        d = ...
            [1	0.3302393550012598	0.0000000000000000
            2	0.1806481606948574	-0.8360311073266358
            3	0.1806481606948574	0.8360311073266358
            4	0.0812743883615744	-0.9681602395076261
            5	0.0812743883615744	0.9681602395076261
            6	0.3123470770400029	-0.3242534234038089
            7	0.3123470770400029	0.3242534234038089
            8	0.2606106964029354	-0.6133714327005904
            9	0.2606106964029354	0.6133714327005904];
    case 10
        d = ...
            [1	0.2955242247147529	-0.1488743389816312
            2	0.2955242247147529	0.1488743389816312
            3	0.2692667193099963	-0.4333953941292472
            4	0.2692667193099963	0.4333953941292472
            5	0.2190863625159820	-0.6794095682990244
            6	0.2190863625159820	0.6794095682990244
            7	0.1494513491505806	-0.8650633666889845
            8	0.1494513491505806	0.8650633666889845
            9	0.0666713443086881	-0.9739065285171717
            10	0.0666713443086881	0.9739065285171717];
end

[~, I] = sort(d(:,3));

W_quad = d(I,2)';
x_quad = d(I,3)';

end