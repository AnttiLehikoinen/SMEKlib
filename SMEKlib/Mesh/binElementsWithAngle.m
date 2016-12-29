function binElements = binElementsWithAngle(elementList, msh, angles, varargin)
%binElements bins elements based on their angular coordinate.
% 
% binElements = binElementsWithAngle(elementList, msh, angles)
% assigns the elements of elementList (referring to the mesh (t, p)) based
% on their mean angular coordinate, to the closest bin on the list
% "angles".
% 
% binElements = binElementsWithAngle(elementList, msh, angles, rs)
% also divides the elements into numel(rs) layers; with rs containing the
% radial coordinates of the inter-layer boundary.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Nbins = numel(angles);

p = msh.p;
t = msh.t;

if numel(varargin) == 1
    rs = varargin{1};
    Nlayers = numel( rs ) + 1;
else
    Nlayers = 1;
    rs = inf;
end

%element center point
Xe = (p(:, t(1,elementList)) + p(:, t(2,elementList)) + p(:, t(3,elementList))) / 3;
elementAngles = atan2(Xe(2,:), Xe(1,:)); elementAngles( elementAngles<0 ) = elementAngles( elementAngles<0 ) + 2*pi;

elementRadii = sqrt( sum(Xe.^2, 1) );

cl = knnsearch(angles', elementAngles')';
binElements = cell(1, Nbins*Nlayers);
for k_angle = 1:Nbins
    for k_layer = 1:Nlayers
        %binElements{1, k_angle} = elementList(cl==k_angle);
        
        if k_layer == 1
            r_low = 0;
        else
            r_low = rs(k_layer-1);
        end
        if k_layer == Nlayers
            r_up = inf;
        else
            r_up = rs(k_layer);
        end
        
        binElements{1, (k_angle-1)*Nlayers + k_layer} = toRow( ...
            elementList( (cl==k_angle) & (elementRadii>r_low) & (elementRadii<r_up) ) );
    end
end

end