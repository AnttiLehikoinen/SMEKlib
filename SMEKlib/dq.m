function xdq = dq(x, angles)

%alpha-beta transformation
xab = 2/3 * [1 -0.5 -0.5;0 sqrt(3)/2 -sqrt(3)/2]*x;

%dq transformation, vectorized
xdq = [cos(angles).*xab(1,:) + sin(angles).*xab(2,:);
    -sin(angles).*xab(1,:) + cos(angles).*xab(2,:)];

end