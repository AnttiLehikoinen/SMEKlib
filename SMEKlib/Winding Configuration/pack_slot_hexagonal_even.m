function [Xc, reff] = pack_slot_hexagonal_even(slotPoints, Ns, rs)
%pack_slot_hexagonal attempts to generate a hexagonal conductor packing.
% 
% [Xc, reff] = pack_slot(slotPoints, Ns, rs)
% tries to generate an approximately uniform packing of Ns round conductors
% of radius rs, inside a slot defined in the 2xN array slotPoints.%
% 
% If required, the effective radius reff is also returned.
% 
% NOTE: rs is assumed to be small enough compared to the slot curvature in
% such a way, that each conductor can only touch two consecutive slot
% boundary segments at once. For slot shapes containing narrow crevices:
% the best of luck.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Xc = [];

%checking orientation
if ispolycw(slotPoints(1,:), slotPoints(2,:))
    slotPoints = fliplr(slotPoints);
end

%getting slot bnd segment normal and tangent vectors
Pstart = slotPoints;
ts = slotPoints(:, [2:end 1]) - slotPoints; %tangent vector
ne = [0 -1;1 0]*ts; ne = bsxfun(@times, ne, 1./sqrt(sum(ts.^2,1))); %unit normal

Nps = size(Pstart, 2);

%coordinate unit vectors
ex = repmat([1;0], 1, Nps); ey = repmat([0;1], 1, Nps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%beginning iterative procedure (division search)
N_iter = 15;
tol = 1e-4;

r_good = rs; r_bad = -1; reff = rs;


for k_iter = 1:N_iter
    %Pstart_eff = Pstart + reff*ne;
    %Pstart_eff = Pstart + reff*(ne + ne(:,[end 1:(end-1)]))/2;
    [Pstart_eff, ts_eff] = internal_contractSlot(Pstart, reff);
    
    didNotFit = 0;

    %finding maximum bottom and top positions for strands
    [Px, ~, ~, k] = internal_findIntersect(0, ex, Pstart_eff, ts_eff);

    inds = find( (k>=(-tol))&(k<=(1+tol)) );
    Xbottom = max(Px(1,inds));
    Xtop = min(Px(1,inds));

    %xcoords = internal_radlinspace(Xbottom, Xtop, reff);
    xcoords = internal_Xtesspace(Xbottom, Xtop, reff);


    Xc_cand = zeros(2, Ns);
    ind_c = 1;
    %clf; hold on;
    for ind_x = 1:numel(xcoords)
        % finding vertical range for strands
        xtemp = repmat([xcoords(ind_x);0], 1, Nps);
        [Py, ~, ~, k] = internal_findIntersect(xtemp, ey, Pstart_eff, ts_eff);
        inds = find( (k>=0)&(k<=1) );
        
        %ycoords = internal_radlinspace(min(Py(2,inds)), max(Py(2,inds)), reff);
        if (ind_x > 1) && any(ycoords)
            ycoords = internal_tesspace(min(Py(2,inds)), max(Py(2,inds)), reff, ycoords(1));
        else
            ycoords = internal_tesspace(min(Py(2,inds)), max(Py(2,inds)), reff, []);
        end

        N_incr = numel(ycoords);
        if (ind_c+N_incr-1) > Ns
            % ran out of strands; quitting
            
            %packing the final layer evenly
            ycoords =aux_intlinspace(min(Py(2,inds)), max(Py(2,inds)), Ns - ind_c + 1);
            
            Xc_cand(1, ind_c:Ns) = xcoords(ind_x);
            Xc_cand(2, ind_c:Ns) = ycoords(1:(Ns-ind_c+1));
            ind_c = Ns+1;
            break;
        else
            Xc_cand(1, ind_c:(ind_c+N_incr-1)) = xcoords(ind_x);
            Xc_cand(2, ind_c:(ind_c+N_incr-1)) = ycoords;
            ind_c = ind_c + N_incr;
        end
        
        %{
        %clf; hold on;
        plot(slotPoints(1,[1:end 1]), slotPoints(2,[1:end 1]), 'b');
        drawConductors(reff, Xc_cand(:,1:(ind_c-1)), 'EdgeColor', 'r'); axis equal;
        plot(xcoords(ind_x)*[1 1], [min(Py(2,inds)) max(Py(2,inds))], 'kx-');
        pause(0.05)
        %}
    end

    if ind_c <= Ns
        didNotFit = 1;
        if k_iter == 1
            error('Conductors do not seem to fit inside the slot!');
        end
    end
    
    % updating bounds
    if didNotFit
        r_bad = reff;
    else
        r_good = reff;
        Xc = Xc_cand;
    end
    
    % plotting for de-bugging purposes
    %{
    figure(6); clf; hold on;
    plotInds = 1:Ns;
    if didNotFit; plotInds = 1:(ind_c-1); end
    
    Np = 20; ap = linspace(0,2*pi,Np)';
    plot(slotPoints(1,[1:end 1]), slotPoints(2,[1:end 1]), 'b');
    XP = bsxfun(@plus, Xc_cand(1,plotInds), reff*cos(ap)); YP = bsxfun(@plus, Xc_cand(2,plotInds), reff*sin(ap));
    plot(XP, YP, 'k'); axis equal;
    
    pause(0.1);
    %}
    
    % updating radius
    if r_bad > 0
        reff = (r_good + r_bad)/2;
    else
        reff = reff * 1.5;
    end
end


end

function xcoords = internal_Xtesspace(Xbottom, Xtop, reff)
%x-coordinates for hexagonal tessellation

if Xbottom > Xtop
    [Xtop, Xbottom] = deal(Xbottom, Xtop);
end

N = floor( (Xtop - Xbottom - reff) / (sqrt(3)*reff) );

xcoords = (Xtop-reff) - (0:(N-1))*(sqrt(3)*reff);


end

function xs = internal_tesspace(y1, y2, reff, yprev)
%y-coordinates for hexagonal tessellation

%{
if y1 > y2
    [y2, y1] = deal(y1, y2);
end
y2 = y2 + reff;
y1 = y1 - reff;
N = floor( (y2-y1)/ (2*reff) );
if mod(N, 2)
    N = N-1;
end
xs = aux_intlinspace(y1, y2, N);
return;
%}


if ~any(yprev)
    %xs = internal_radlinspace(x1, x2, reff);
    yprev = -reff;
    %return;
end
%yprev = -reff;

if y1 > y2
    [y2, y1] = deal(y1, y2);
end

N_up = floor( abs(y2 - yprev - reff) / (2*reff) );
%N_down = floor( abs(y1 - yprev - reff) / (2*reff) );
N_down = floor( abs(yprev+reff - y1) / (2*reff) );

xs = yprev + ((-N_down):N_up)*reff*2 + reff;

TOL = 0.2e-3;
xs = xs( abs(xs) > TOL );

end


function [Pstart_eff, ts_eff] = internal_contractSlot(Pstart, reff)
% tries to move the slot boundaries in by the amount reff
% algorithm used:
% http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

ts = Pstart(:, [2:end 1]) - Pstart; %tangent vectors
ne = [0 -1;1 0]*ts; ne = bsxfun(@times, ne, 1./sqrt(sum(ts.^2,1))); %unit inner normal vectors

%boundary segments translated inwards by rs
P1 = Pstart + reff*ne; t1 = ts;
P2 = P1(:, [end 1:(end-1)]); t2 = ts(:, [end 1:(end-1)]);

    %nested function for 2D cross-product
    function o = aux_cross(V1, V2)
        o = V1(1,:).*V2(2,:) - V2(1,:).*V1(2,:);
    end

%finding intersection points of
% type P1 + c1*t1 = P2 + c2*t2
denom = aux_cross(t1, t2);
c1 = aux_cross(P2-P1, t2) ./ denom;


Pstart_eff = P1 + bsxfun(@times, t1, c1);
Pstart_eff(:, abs(denom)<1e-6) = P1(:, abs(denom)<1e-6);

ts_eff = Pstart_eff(:, [2:end 1]) - Pstart_eff;


end

function [p1, k1, p2, k2] = internal_findIntersect(P1, t1, P2, t2)
%finds the intersects of two sets of lines

%figure(5); plot(P2(1,:), P2(2,:), 'k.');

M = [dotProduct(t1, t1); dotProduct(t1, t2); dotProduct(t2, t2)];
detM = -M(1,:).*M(3,:) + M(2,:).^2;

F = [dotProduct(P2-P1, t1); dotProduct(P2-P1,t2)];

coeffs = mappingTimesVector(F, 1, 0, [M(1:2,:); -M(2:3,:)], [], detM);
k1 = coeffs(1,:);
k2 = coeffs(2,:);

p1 = P1 + bsxfun(@times, k1, t1);
p2 = P2 + bsxfun(@times, k2, t2);

end

function xs = internal_radlinspace(x1, x2, reff)

lint = x2 - x1; %interval length
Nc = floor( abs(lint) / (2*reff) ); %number of strands fitting
lfree = (abs(lint) - 2*Nc*reff) / (Nc+1);

if Nc == 1
    xs = (x1+x2)/2;
else
    xs = linspace(x1 + sign(lint)*lfree, x2 - sign(lint)*lfree, Nc);
end

end