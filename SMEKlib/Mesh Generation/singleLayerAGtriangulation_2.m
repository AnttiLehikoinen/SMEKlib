function tag = singleLayerAGtriangulation_2(msh, n_ag_s, n_ag_r)
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

tag = zeros(3, 1000); ne = 0;
function addt(t)
    % nested function for incrementing the array
    if (ne+size(t,2)) > size(tag,2)
        nadd = min(2*size(tag,2), size(t,2));
        tag = [tag zeros(2,nadd)];
    end
    tag(:,(ne+1):(ne+size(t,2))) = t;
    ne = ne + size(t,2);
end

%triangulating
rprev = 1; N_ag_s = numel(n_ag_s);
sprev = 1; N_ag_r = numel(n_ag_r);
while true
    %triangulation for the next 4 air-gap nodes (2 candidate triangles)
    TR = triangulation([1 2 3; 1 3 4], [msh.p(:, n_ag_s(sprev:(sprev+1))) msh.p(:,n_ag_r(rprev:(rprev+1)))]');
    
    %computing the radii of in- (rin) and circumcircles (rout)
    [~,rin1] = incenter(TR, 1); [~,rout1] = circumcenter(TR, 1);
    [~,rin2] = incenter(TR, 2); [~,rout2] = circumcenter(TR, 2);
    
    %pt = [msh.p(:, n_ag_s(sprev:(sprev+1))) msh.p(:,n_ag_r(rprev:(rprev+1)))];
    %[rout1, rout2] = circumcircles(pt);
    
    %choosing the triangle with the better condition indicator (rout/rin)
    if rout1/rin1 < rout2/rin2
        %tag(:,k) = [n_ag_s(sprev:(sprev+1)) n_ag_r(rprev)]';
        addt( [n_ag_s(sprev:(sprev+1)) n_ag_r(rprev)]' );
        sprev = sprev+1;
    else
        %tag(:,k) = [n_ag_s(sprev) n_ag_r(rprev:(rprev+1))]';
        addt( [n_ag_s(sprev) n_ag_r(rprev:(rprev+1))]' );
        rprev = rprev+1;
    end
    
    %checking if running out of nodes on either boundary
    if sprev == N_ag_s
        while rprev < N_ag_r
            addt( [n_ag_s(sprev) n_ag_r(rprev:(rprev+1))]' );
            rprev = rprev+1;
        end
        break;
    end
    if rprev == N_ag_r
        while sprev < N_ag_s
            addt( [n_ag_s(sprev:(sprev+1)) n_ag_r(rprev)]' );
            sprev = sprev+1;
        end
        break;
    end
end
tag = tag(:,1:ne);

end

function [rout1, rout2] = circumcircles(pt)
a = norm(pt(:,2) - pt(:,1)); 
b = norm(pt(:,3)-pt(:,2)); 
c = norm(pt(:,1)-pt(:,3));
a2 = norm(pt(:,4)-pt(:,3));
b2 = norm(pt(:,4)-pt(:,1));

A1 = abs(det([pt(:,2)-pt(:,1) pt(:,3)-pt(:,1)]));
A2 = abs(det([pt(:,3)-pt(:,1) pt(:,4)-pt(:,1)]));
rout1 = a*b*c/(2*A1);
rout2 = a2*b2*c/(2*A2);
end

