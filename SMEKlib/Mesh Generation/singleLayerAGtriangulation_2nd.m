function [tag, pextra] = singleLayerAGtriangulation_2nd(msh, n_ag_s, n_ag_r)
%singleLayerAGtriangulation Single-layer air-gap triangulation.
%
% singleLayerAGtriangulation(msh, n_ag_s, n_ag_r) tries to generate a
% single-layer triangulation between the stator (n_ag_s) and rotor (n_ag_r)
% air-gap boundary nodes.
%
% (c) 2017 Antti Lehikoinen / Aalto University

%sorting ag nodes
%[~,I] = sort( atan2(msh.p(2,n_ag_s), msh.p(1,n_ag_s)) ); n_ag_s = n_ag_s(I);
%[~,I] = sort( atan2(msh.p(2,n_ag_r), msh.p(1,n_ag_r)) ); n_ag_r = n_ag_r(I);


tag = zeros(6, 1000); ne = 0;
    function addt(t)
        % nested function for incrementing the array
        if (ne+size(t,2)) > size(tag,2)
            nadd = min(2*size(tag,2), size(t,2));
            tag = [tag zeros(6,nadd)];
        end
        tag(:,(ne+1):(ne+size(t,2))) = t;
        ne = ne + size(t,2);
    end

pextra = zeros(2, 1000); np = 0;
    function addp(p)
        % nested function for incrementing the array
        if (np+size(p,2)) > size(pextra,2)
            nadd = min(2*size(pextra,2), size(p,2));
            pextra = [pextra zeros(2,nadd)];
        end
        pextra(:,(np+1):(np+size(p,2))) = p;
        np = np + size(p,2);
    end

%triangulating
Np = size(msh.p,2);
rprev = 1; N_ag_s = numel(n_ag_s);
sprev = 1; N_ag_r = numel(n_ag_r);

figure(10); clf; hold on;
while true
    %triangulation for the next 4 air-gap nodes (2 candidate triangles)
    %TR = triangulation([1 2 3; 1 3 4], [msh.p(:, n_ag_s(sprev:(sprev+1))) msh.p(:,n_ag_r(rprev:(rprev+1)))]');
    
    %computing the radii of in- (rin) and circumcircles (rout)
    %[~,rin1] = incenter(TR, 1); [~,rout1] = circumcenter(TR, 1);
    %[~,rin2] = incenter(TR, 2); [~,rout2] = circumcenter(TR, 2);
    
    pt = [msh.p(:, n_ag_s(sprev:2:(sprev+2))) msh.p(:,n_ag_r(rprev:2:(rprev+2)))];
    [rout1, rout2, rin1, rin2] = circumcircles(pt);
    
    %choosing the triangle with the better condition indicator (rout/rin)
    if rout1/rin1 < rout2/rin2
        %tag(:,k) = [n_ag_s(sprev:(sprev+1)) n_ag_r(rprev)]';
        tt = [n_ag_s(sprev:2:(sprev+2)) n_ag_r(rprev) ...
            n_ag_s(sprev+1) ...
            Np+np+2 ...
            Np+np+1]';
        addt( tt );
        
        pt = 0.5*(msh.p(:, n_ag_s(sprev))+msh.p(:,n_ag_r(rprev)));
        
        addp(pt);
        
        %triplot(tt(1:3)', msh.p(1,:), msh.p(2,:)); drawnow;
        
        sprev = sprev+2;
    else
        %tag(:,k) = [n_ag_s(sprev) n_ag_r(rprev:(rprev+1))]';
        %addt( [n_ag_s(sprev) n_ag_r(rprev:(rprev+1))]' );
        tt = [n_ag_r(rprev:2:(rprev+2)) n_ag_s(sprev) ...
            n_ag_r(rprev+1) ...
            Np+np+2 ...
            Np+np+1]';
        addt( tt );

        pt = 0.5*(msh.p(:, n_ag_r(rprev))+msh.p(:,n_ag_s(sprev)));
        
        addp(pt);
        rprev = rprev+2;
    end
    
    %checking if running out of nodes on either boundary
    if sprev >= (N_ag_s-1)
        while rprev < N_ag_r
            tt = [n_ag_r(rprev:2:(rprev+2)) n_ag_s(sprev) ...
            n_ag_r(rprev+1) ...
            Np+np+2 ...
            Np+np+1]';
        addt( tt );
        
        addp( 0.5*(msh.p(:, n_ag_r(rprev))+msh.p(:,n_ag_s(sprev))) );
        rprev = rprev+2;
        end
        break;
    end
    if rprev >= (N_ag_r-1)
        while sprev < N_ag_s
            %addt( [n_ag_s(sprev:(sprev+1)) n_ag_r(rprev)]' );
            tt = [n_ag_s(sprev:2:(sprev+2)) n_ag_r(rprev) ...
                n_ag_s(sprev+1) ...
                Np+np+2 ...
                Np+np+1]';
            addt( tt );
            
            addp( 0.5*(msh.p(:, n_ag_s(sprev))+msh.p(:,n_ag_r(rprev))) );
            
            %triplot(tt(1:3)', msh.p(1,:), msh.p(2,:)); drawnow;
            
            sprev = sprev+2;
        end
        break;
    end
end
addp( 0.5*(msh.p(:, n_ag_s(end))+msh.p(:,n_ag_r(end))) );

tag = tag(:,1:ne);
pextra = pextra(:,1:np);

end

function [rout1, rout2, rin1, rin2] = circumcircles(pt)
a = norm(pt(:,2) - pt(:,1));
b = norm(pt(:,3)-pt(:,2));
c = norm(pt(:,1)-pt(:,3));
a2 = norm(pt(:,4)-pt(:,3));
b2 = norm(pt(:,4)-pt(:,1));

A1 = abs(det([pt(:,2)-pt(:,1) pt(:,3)-pt(:,1)]));
A2 = abs(det([pt(:,3)-pt(:,1) pt(:,4)-pt(:,1)]));
rout1 = a*b*c/(2*A1);
rout2 = a2*b2*c/(2*A2);

s1 = (a+b+c)/2; s2 = (a2+b2+c)/2;
rin1 = A1/s1; rin2 = A2/s2;

end

