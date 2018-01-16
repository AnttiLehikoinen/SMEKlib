function coeffs = preCompute_referenceShapeFunctions(elementType, shapefunction)

if elementType==Elements.prism && shapefunction==ShapeFunctions.nedelec
    %setting up a prismatic mesh of one reference element
    m2 = MeshBase([0 0;1 0;0 1]', [1 2 3]');
    m3 = ExtrudedPrismMesh(m2, [0 1]);
    
    %edge tangents
    Tau = m3.nodes(:, m3.edges(2,:)) - m3.nodes(:, m3.edges(1,:)); %global
    %shifting to element-wise ordering
    temp = m3.elements2edges(:,1)';
    Tau = bsxfun(@times, Tau(:, abs(temp)), sign(temp) );
    
    %starting points of edges
    Tau0 = m3.nodes(:, m3.edges(1, abs(temp)));
    Tau0(:, temp<0) = m3.nodes(:, m3.edges(2, -temp(temp<0)));
    
    Pcoeffs = get_NedelecPolynomialCoefficients(elementType, shapefunction);
    bt = PolynomialBasis3D.getBasisType(elementType, shapefunction);
    C = zeros(9);
    for kedge = 1:9
        for kp = 1:9
            C(kedge, kp) = computeEdgeTangentIntegral(Pcoeffs{kp}, Tau(:,kedge), ...
                Tau0(:,kedge), bt);
        end
    end
    coeffs_b = transpose( C\eye(9) );
    coeffs = cell(1, 9);
    for k = 1:9
        coeffs{k} = zeros(3, 6);
        for kb = 1:9
            coeffs{k} = coeffs{k} + coeffs_b(k, kb) * Pcoeffs{kb};
        end
    end
else
    error('Not yet implemented.');
end

end

function I = computeEdgeTangentIntegral(coeffs, tau, tau0, bt)

[x, w] = get_1DQuadPoints(3); %correct??
x = 0.5 + 0.5*x;
w = 0.5*w; %shifting to integral over [0,1]

I = 0;
tau_n = tau / norm(tau);
for k = 1:numel(w)
    X = tau0 + x(k)*tau;
    I = I + w(k) * tau_n' * (coeffs * PolynomialBasis3D.bteval(X, bt));
end

end