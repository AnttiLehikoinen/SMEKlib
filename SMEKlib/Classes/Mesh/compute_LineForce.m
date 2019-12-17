function F = compute_LineForce(Br, msh1_top, msh1_bottom, msh2_top, msh2_bottom)

%quadrature points
[xq, wq] = msh1_top.quadpoints(1);

%mappings
[F1t, f1t] = msh1_top.getMappingMatrix(); det1t = matrixDeterminant(F1t);
[F1b, f1b] = msh1_bottom.getMappingMatrix(); det1b = matrixDeterminant(F1b);
[F2t, f2t] = msh2_top.getMappingMatrix(); det2t = matrixDeterminant(F2t);
[F2b, f2b] = msh2_bottom.getMappingMatrix(); det2b = matrixDeterminant(F2b);

%looping
F = zeros(3, 1);
mu0 = pi*4e-7;

for kq = 1:numel(wq)
    %global points
    x1t = matrixTimesVector(F1t, xq(:,kq), false, false) + f1t;
    x1b = matrixTimesVector(F1b, xq(:,kq), false, false) + f1b;
    
    
    Hs_top = zeros(3, size(msh1_top.t, 2));
    Hs_bottom = zeros(3, size(msh1_top.t, 2));
    
    for kq2 = 1:numel(wq)
        x2t = matrixTimesVector(F2t, xq(:,kq2), false, false) + f2t;
        x2b = matrixTimesVector(F2b, xq(:,kq2), false, false) + f2b;
    
        for ke2 = 1:size(msh2_top.t, 2)
            dx = x1t - x2t(:, ke2);
            Hs_top = Hs_top + dx ./ sum(dx.^2, 1) * det2t(ke2)*wq(kq2);

            dx = x1t - x2b(:, ke2); 
            Hs_top = Hs_top - dx ./ sum(dx.^2, 1) * det2b(ke2)*wq(kq2);

            dx = x1b - x2t(:, ke2); 
            Hs_bottom = Hs_bottom - dx ./ sum(dx.^2, 1) * det2t(ke2)*wq(kq2);

            dx = x1b - x2b(:, ke2); 
            Hs_bottom = Hs_bottom + dx ./ sum(dx.^2, 1) * det2b(ke2)*wq(kq2);
        end
    end
    
    %H = wq(1) * H / (4*pi) * Br;
    %F = F + wq(1)*H*Br*detF(k1) / mu0;
    
    F = F + Br^2*wq(kq) * ( sum(Hs_top / (4*pi) .* det1t, 2) + ...
        sum(Hs_bottom / (4*pi) .* det1b, 2) );
end

F = 2*F / mu0; %2-coefficient excluded from inner integration

end
