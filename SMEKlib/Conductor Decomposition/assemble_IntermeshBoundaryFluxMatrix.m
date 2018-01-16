function S = assemble_IntermeshBoundaryFluxMatrix(msh_master, msh_slave, boundaryData, v_master, S)
% assume two meshes "master" and "slave", joined together at a boundary
% defined in boundaryData.
% This function returns a matrix S s.t. the vector
% B_m = S * A_s
% contains the entries
% [B_m]_i = Int{ Phi_m * dA_s/dN }
% where Phi_m are the shape functions of the master mesh, while A_s is the
% discretized vector potential solution in the slave mesh

% boundaryData: a 3 x Nb matrix describing the boundary between the meshed
% in the following segmented way
%   - boundaryData(1,:) = master element this segment belongs to
%   - boundaryData(2:3,:) = start node coordinates in the master mesh
%   - boundaryData(4:5,:) = end node coordinates
%   - boundaryData(6:10,:) = same data for slave mesh

N_segments = size(boundaryData, 2);

%determining mappings
[F_master, F0_master] = mappingTerms(msh_master, boundaryData(1,:)); detF_master = mappingDeterminant(F_master);
[F_slave, F0_slave] = mappingTerms(msh_slave, boundaryData(6,:)); detF_slave = mappingDeterminant(F_slave);

%getting reference 1D quad points
N_quad = 5;
[x_quad, W_quad] = internal_getQuadPoints(N_quad);
W_quad = W_quad/2; %scaling for interval length

%getting segment start and end points in respective reference elements
Xref_start_master = mappingTimesVector(boundaryData(2:3,:), 1, 0, F_master, -F0_master, detF_master);
Xref_end_master = mappingTimesVector(boundaryData(4:5,:), 1, 0, F_master, -F0_master, detF_master);

Xref_start_slave = mappingTimesVector(boundaryData(7:8,:), 1, 0, F_slave, -F0_slave, detF_slave);
Xref_end_slave = mappingTimesVector(boundaryData(9:10,:), 1, 0, F_slave, -F0_slave, detF_slave);

%getting slave segment normal vectors
ne_slave = [0 1;-1 0] * (boundaryData(9:10,:) - boundaryData(7:8,:));
le_slave = sqrt( sum(ne_slave.^2, 1) );
ne_slave = bsxfun(@times, ne_slave, 1./le_slave);
le_master = sqrt(sum((boundaryData(4:5,:) - boundaryData(2:3,:)).^2,1));

%defining function handles and indices
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] };
fh_master = @(ind_n, X)( phi_ref{ind_n}*X );
fh_slave = @(ind_n, X)( mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F_slave, [], detF_slave) );
NO_f1 = 3; NO_f2 = 3;
INDM_1 = msh_master.t(:, boundaryData(1,:));
INDM_2 = msh_slave.t(:, boundaryData(6,:));

%initializing quadrature calculation
E = cell(NO_f1, NO_f2); [E{:}] = deal( zeros(1, N_segments) );

%calculating entries
for k_quad = 1:N_quad
    X_quad_master = [(Xref_end_master - Xref_start_master)*x_quad(k_quad) / 2 + ...
        (Xref_end_master + Xref_start_master) / 2;
        ones(1, N_segments)];
    X_quad_slave = [(Xref_end_slave - Xref_start_slave)*x_quad(k_quad) / 2 + ...
        (Xref_end_slave + Xref_start_slave) / 2;
        ones(1, N_segments)];
    for ind_1 = 1:NO_f1
        for ind_2 = 1:NO_f2            
            E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) * ...
                fh_master(ind_1, X_quad_master) .* dotProduct( fh_slave(ind_2, X_quad_slave), ne_slave );
            %E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) * fh_master(ind_1, X_quad_master);
            %E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) .* dotProduct( fh_slave(ind_2, X_quad_slave), ne_slave );
        end
    end
end

%adding entries to matrix struct
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        S = sparseAdd( INDM_1(ind_1,:), INDM_2(ind_2,:), ...
            v_master .* E{ind_1,ind_2}.*le_slave, S);
    end
end

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