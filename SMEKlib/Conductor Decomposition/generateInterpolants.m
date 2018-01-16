function [interpolants, Pclosest] = generateInterpolants( PM, PS, bnd_edges_master, slavePoints, PsFun, varargin )
Pclosest = PS(:, slavePoints);
if isa(bnd_edges_master, 'cell')
    %[interpolants, Pclosest] = internal_varyingOrderInt( PM, PS, bnd_edges_master, slavePoints, PsFun );
    [interpolants, Pclosest] = internal_varyingOrderInt_2( PM, PS, bnd_edges_master, slavePoints, PsFun );
    return;
end

if ~isempty(varargin)
    Norder = varargin{1};
    interpolants = internal_highOrderInt( PM, PS, bnd_edges_master, slavePoints, PsFun, Norder );
    return;
end

if size(bnd_edges_master,1) == 2
    interpolants = internal_firstOrder( PM, PS, bnd_edges_master, slavePoints, PsFun );
elseif size(bnd_edges_master,1) == 3
    interpolants = internal_secondOrder( PM, PS, bnd_edges_master, slavePoints, PsFun );
end


end

function interpolants = internal_firstOrder( PM, PS, bnd_edges_master, slavePoints, PsFun )

P_slave = PsFun( PS(:, slavePoints) );

%closest master edge to each slave point
closest_s2m = closestSegmentToNode( P_slave, ...
   PM(:, bnd_edges_master(1,:)),PM(:, bnd_edges_master(2,:)) );

Npoints = 2;
Np_slave = size(closest_s2m, 2);

I = zeros(Np_slave*Npoints, 1); J = I; E = I;
ri = 1;
for kp_slave = 1:Np_slave
    me = closest_s2m(1, kp_slave);
    m_start = PM(:, bnd_edges_master(1,me));
    m_end = PM(:, bnd_edges_master(2,me));
    l_edge = norm( m_start - m_end );
    l_point = norm( closest_s2m(2:3, kp_slave) - m_start );
    alpha = l_point / l_edge;
    
    I(ri:(ri+Npoints-1)) = slavePoints(kp_slave) * ones(Npoints, 1);
    J(ri:(ri+Npoints-1)) = [bnd_edges_master(1,me); bnd_edges_master(2,me)];
    E(ri:(ri+Npoints-1)) = [1-alpha alpha];
    
    ri = ri + Npoints;
end

interpolants = [I J E];

end

function interpolants = internal_secondOrder( PM, PS, bnd_edges_master, slavePoints, PsFun )

N_newton = 5;

%reference shape functions
phi_ref = {[1 -3 2]; [0 -1 2]; [0 4 -4] };

P_slave = PsFun( PS(:, slavePoints) );

%%%%%%%%%%%%%%
%testing
%phi_ref{:}
%C_ref = internal_getRefShapeFunctions(2)
%tb_fun = @(t, dN)( internal_evaluatePolyBasis(t, 2, dN) );
%Pe_fun = @(tb)( internal_evaluateEdgeCoordinates(tb, C_ref, PM,bnd_edges_master) );

%linear edge approximation
t_master = PM(:, bnd_edges_master(2,:)) - PM(:, bnd_edges_master(1,:));
tprod = dotProduct(t_master, t_master);
P1 = PM(:, bnd_edges_master(1,:));

%edge coordinate function
tbfun = @(t)( [t.^0;t.^1;t.^2] );
d_tbfun = @(t)( [0*t;t.^0;2*t.^1] );
dd_tbfun = @(t)( [0*t;0*t;2*t.^0] );
fun_Pe = @(tb)( bsxfun(@times, PM(:, bnd_edges_master(1,:)), (phi_ref{1}*tb)) + ...
    bsxfun(@times, PM(:, bnd_edges_master(2,:)), (phi_ref{2}*tb)) + ...
    bsxfun(@times, PM(:, bnd_edges_master(3,:)), (phi_ref{3}*tb)) );

Npoints = 3;
Np_slave = size(P_slave, 2);
Ne_master = size(bnd_edges_master, 2);

I = zeros(Np_slave*Npoints, 1); J = I; E = I;
ri = 1;
%going through slave points
for kp_slave = 1:Np_slave
    %getting starting point for each segment
    Ps = repmat(P_slave(:, kp_slave), 1, Ne_master);
    
    tvalues = ( inner_dotProduct(Ps, t_master) - dotProduct(P1, t_master) ) ./ tprod;
    tvalues = max(0, min(tvalues, 1));
    
    %running a few Newton iterations to improve guess
    for k_Newton = 1:N_newton
        Pe = fun_Pe( tbfun(tvalues) );
        dPe = fun_Pe( d_tbfun(tvalues) );
        ddPe = fun_Pe( dd_tbfun(tvalues) );
        
        dF = 2*dotProduct(dPe, Pe-Ps);
        ddF = 2*dotProduct(ddPe, Pe-Ps) + 2*dotProduct(dPe, dPe);
        
        tvalues = tvalues - dF./ddF;
        tvalues = max(0, min(tvalues, 1));
    end
    
    Pe = fun_Pe( tbfun(tvalues) );
    dists = dotProduct(Pe,Pe - 2*Ps);
    [~,me] = min(dists);
    tb = tbfun( tvalues(me) );
    
    %tb2 = tb_fun(tvalues(:,me), 0);     disp([tb tb2])
    %Pe2 = Pe_fun( tb ); disp([Pe(:,me) Pe2(:,me)])    
    %disp( dotProduct(Pe(:,me),Pe(:,me) - 2*Ps(:,me)) + dotProduct(Ps(:,me), Ps(:,me)) )
    
    I(ri:(ri+Npoints-1)) = slavePoints(kp_slave) * ones(Npoints, 1);
    J(ri:(ri+Npoints-1)) = bnd_edges_master(:, me);
    E(ri:(ri+Npoints-1)) = [phi_ref{1}*tb; phi_ref{2}*tb; phi_ref{3}*tb];
    
    ri = ri + Npoints;
end

interpolants = [I J E];

end

function [interpolants, P_closest] = internal_varyingOrderInt( PM, PS, EdgeSets, slavePoints, PsFun )

N_newton = 5;

Nsets = numel( EdgeSets );

%generating stuff for each edge
t_lin = cell(Nsets, 1); %linear estimate of edge tangent vectors
t_prod = cell(Nsets, 1);
P_start = cell(Nsets, 1);

N_edgeOrder = zeros(1, Nsets);
C_ref = cell(Nsets, 1);
tbfun = cell(Nsets, 1);
Pe_fun = cell(Nsets, 1);

for ke = 1:Nsets
    N_edgeOrder(ke) = size(EdgeSets{ke}, 1) - 1;
    
    t_lin{ke} = PM(:, EdgeSets{ke}(N_edgeOrder(ke)+1,:)) - PM(:, EdgeSets{ke}(1,:));
    t_prod{ke} = dotProduct(t_lin{ke}, t_lin{ke});
    P_start{ke} = PM(:, EdgeSets{ke}(1,:));
    
    C_ref{ke} = internal_getRefShapeFunctions(N_edgeOrder(ke));
    
    tbfun{ke} = @(t, dN)( internal_evaluatePolyBasis(t, N_edgeOrder(ke), dN) );
    Pe_fun{ke} = @(tb)( internal_evaluateEdgeCoordinates(tb, C_ref{ke}, PM, EdgeSets{ke}) );
end

%moving slave points if required
P_slave = PsFun(PS);
P_closest = zeros(2, numel(slavePoints));

Npoints = max(N_edgeOrder) + 1;
Np_slave = numel(slavePoints);

I = zeros(Np_slave*Npoints, 1); J = I; E = I;
ri = 1;

%going through slave points
for kp_slave = 1:Np_slave
    
    %going through edges of different order
    minDists = zeros(Nsets, 1); minInds = minDists; min_tb = cell(Nsets, 1); min_Pe = zeros(2, Nsets);
    for ke = 1:Nsets
        %getting starting point for each segment
        Ps = repmat(P_slave(:, slavePoints(kp_slave)), 1, size(EdgeSets{ke}, 2));

        tvalues = ( inner_dotProduct(Ps, t_lin{ke}) - dotProduct(P_start{ke},t_lin{ke}) ) ./ t_prod{ke};
        tvalues = max(0, min(tvalues, 1));

        %running a few Newton iterations to improve guess
        for k_Newton = 1:N_newton
            Pe = Pe_fun{ke}( tbfun{ke}(tvalues, 0) );
            dPe = Pe_fun{ke}( tbfun{ke}(tvalues, 1) );
            ddPe = Pe_fun{ke}( tbfun{ke}(tvalues, 2) );

            dF = 2*dotProduct(dPe, Pe-Ps);
            ddF = 2*dotProduct(ddPe, Pe-Ps) + 2*dotProduct(dPe, dPe);

            tvalues = tvalues - dF./ddF;
            tvalues = max(0, min(tvalues, 1));
        end
        tb_final = tbfun{ke}(tvalues, 0);
        
        Pe = Pe_fun{ke}( tb_final );
        dists = dotProduct(Pe,Pe - 2*Ps);
        [minD, minE] = min(dists);
        
        %saving best point on this edge set
        minDists(ke) = minD; minInds(ke) = minE; min_tb{ke} = tb_final(:, minE); min_Pe(:,ke) = Pe(:, minE);
    end
    
    %finding globally best approximation
    [~,ind_set] = min(minDists);
    N_incr = N_edgeOrder(ind_set) + 1;
    
    %saving closest point
    P_closest(:, kp_slave) = min_Pe(:, ind_set);
    
    %adding entries
    I(ri:(ri+N_incr-1)) = slavePoints(kp_slave) * ones(N_incr, 1);    
    J(ri:(ri+N_incr-1)) = EdgeSets{ind_set}(:, minInds(ind_set));
    E(ri:(ri+N_incr-1)) = C_ref{ind_set} * min_tb{ind_set};
    
    ri = ri + N_incr;
    
end

I = I(1:(ri-1)); J = J(1:(ri-1)); E = E(1:(ri-1));

interpolants = [I J E];

end

function [interpolants, P_closest] = internal_varyingOrderInt_2( PM, PS, EdgeSets, slavePoints, PsFun )
%now with a hopefully faster implementation

%some quirks for this method
if iscolumn(slavePoints)
    slavePoints = slavePoints';
end
    

N_iter = 7;

Nsets = numel( EdgeSets );

%generating stuff for each edge
t_lin = cell(Nsets, 1); %linear estimate of edge tangent vectors
t_prod = cell(Nsets, 1);
P_start = cell(Nsets, 1);

N_edgeOrder = zeros(1, Nsets);
C_ref = cell(Nsets, 1);
tbfun = cell(Nsets, 1);
Pe_fun = cell(Nsets, 1);

for ke = 1:Nsets
    N_edgeOrder(ke) = size(EdgeSets{ke}, 1) - 1;
    
    t_lin{ke} = PM(:, EdgeSets{ke}(N_edgeOrder(ke)+1,:)) - PM(:, EdgeSets{ke}(1,:));
    t_prod{ke} = dotProduct(t_lin{ke}, t_lin{ke});
    P_start{ke} = PM(:, EdgeSets{ke}(1,:));
    
    C_ref{ke} = internal_getRefShapeFunctions(N_edgeOrder(ke));
    
    tbfun{ke} = @(t, dN)( internal_evaluatePolyBasis(t, N_edgeOrder(ke), dN) );
    %Pe_fun{ke} = @(tb)( internal_evaluateEdgeCoordinates(tb, C_ref{ke}, PM, EdgeSets{ke}) );
end

%moving slave points if required
P_slave = PsFun(PS(:, slavePoints));
P_closest = zeros(2, numel(slavePoints));

Npoints = max(N_edgeOrder) + 1;
Np_slave = numel(slavePoints);

%going through all edges first
I_initial = zeros(Npoints, Np_slave); J_initial = I_initial; E_initial = I_initial;
minDists = 1e16*ones(1, Np_slave);

for kset = 1:Nsets
    for ke = 1:size(EdgeSets{kset},2)
        %closest point on this edge (linear approximation) to each slave
        %point
        tvalues = ( internal_dotProduct(t_lin{kset}(:,ke), P_slave) - internal_dotProduct(t_lin{kset}(:,ke), P_start{kset}(:,ke)) ) / t_prod{kset}(ke);
        tvalues = max(0, min(tvalues, 1));
        
        %tangent method
        for k_iter = 1:N_iter*0
            %tangent vectors at the estimated points
            T = internal_evaluateEdgeCoordinates(tbfun(tvalues,1), C_ref{kset}, PM, EdgeSets{kset}(:,ke));
            Tnorm = sqrt( dotProduct(T,T) );
            T = bsxfun(@times, T, 1./Tnorm);
            
            %estimated closest point
            Pest = internal_evaluateEdgeCoordinates(tbfun(tvalues,0), C_ref{kset}, PM, EdgeSets{kset}(:,ke));
            
            %new estimated closest point on the tangent vector
            ttemp = (dotProduct(T, P_slave) - dotProduct(T,Pest)) ./ Tnorm;
            
            tvalues = tvalues + ttemp;
            tvalues = max(0, min(tvalues, 1));
        end
        
        %Newton's method
        for k_iter = 1:N_iter            
            Pe = internal_evaluateEdgeCoordinates2(tbfun{kset}(tvalues, 0), C_ref{kset}, PM, EdgeSets{kset}(:,ke));
            dPe = internal_evaluateEdgeCoordinates2(tbfun{kset}(tvalues, 1), C_ref{kset}, PM, EdgeSets{kset}(:,ke));
            ddPe = internal_evaluateEdgeCoordinates2(tbfun{kset}(tvalues, 2), C_ref{kset}, PM, EdgeSets{kset}(:,ke));

            dF = 2*dotProduct(dPe, Pe-P_slave);
            ddF = 2*dotProduct(ddPe, Pe-P_slave) + 2*dotProduct(dPe, dPe);

            tvalues = tvalues - dF./ddF;
            tvalues = max(0, min(tvalues, 1));
        end
        
        %final approximated point
        tb_final = tbfun{kset}(tvalues, 0);        
        Pe = internal_evaluateEdgeCoordinates2(tb_final, C_ref{kset}, PM, EdgeSets{kset}(:,ke));
        dists = dotProduct(Pe,Pe - 2*P_slave);
        
        %finding where we improved
        inds_better = find( dists < minDists );
        minDists(inds_better) = dists( inds_better );
        
        %saving improved results
        if any(inds_better)
            sNodes = slavePoints( inds_better ); 
            mNodes = EdgeSets{kset}(:,ke); 
            coeffs = C_ref{kset}*tb_final(:,inds_better);
            I_initial(:,inds_better) = 0*I_initial(:,inds_better);
            
            I_initial(1:(N_edgeOrder(kset)+1), inds_better) = repmat(sNodes, N_edgeOrder(kset)+1, 1);
            J_initial(1:(N_edgeOrder(kset)+1), inds_better) = repmat(mNodes, 1, numel(inds_better));
            E_initial(1:(N_edgeOrder(kset)+1), inds_better) = coeffs;

            P_closest(:,inds_better) = Pe(:, inds_better);
        end
    end
end

%final construction
I = zeros(Np_slave*Npoints, 1); J = I; E = I;
ri = 1;

for kp_slave = 1:Np_slave
    N_incr = numel(find(I_initial(:,kp_slave)));
    I(ri:(ri+N_incr-1)) = I_initial(1:N_incr, kp_slave);  
    J(ri:(ri+N_incr-1)) = J_initial(1:N_incr, kp_slave);  
    E(ri:(ri+N_incr-1)) = E_initial(1:N_incr, kp_slave);  
    
    ri = ri + N_incr;
end

I = I(1:(ri-1)); J = J(1:(ri-1)); E = E(1:(ri-1));

interpolants = [I J E];

end

function E = internal_dotProduct(x, Y)
% dot product between one vector and several vectors
E = x(1)*Y(1,:) + x(2,:)*Y(2,:);
end

function C = internal_getRefShapeFunctions(Norder)
%returns coefficients of the Lagrance shape functions of order Norder,
%evaluated at the reference edge

tpoints = linspace(0, 1, Norder + 1);

M = bsxfun(@power, tpoints', 0:Norder);


C = transpose( M \ eye(Norder+1) );

%valueMatrix = eye(Norder+1);
%C = transpose( M \ valueMatrix([1 Norder+1 2:Norder],:) );

end

function tb = internal_evaluatePolyBasis(t, Norder, derivOrder)
%returns the univariate polynomial basis of specified order and derivative
%order

if derivOrder == 0
    tb = bsxfun(@power, t, (0:Norder)' );
elseif derivOrder == 1
    %tb = [diag(1:Norder) * bsxfun(@power, t, (0:(Norder-1))' ); 
    %    zeros(1, numel(t))];
    tb = [zeros(1, numel(t));
        diag(1:Norder) * bsxfun(@power, t, (0:(Norder-1))' )];
else
    tb = [zeros(2, numel(t)); 
        diag(1:(Norder-1))*diag(2:Norder) * bsxfun(@power, t, (0:(Norder-2))' )];
end


end

function Pe = internal_evaluateEdgeCoordinates(tb, C, P, edges)

Pe = zeros(2, size(edges,2));

for ke = 1:size(edges,1)
    %size( bsxfun(@times, P(:, edges(ke,:)), C(ke,:)*tb) )
    Pe = Pe + bsxfun(@times, P(:, edges(ke,:)), C(ke,:)*tb);
end

end

function Pe = internal_evaluateEdgeCoordinates2(tb, C, P, edge)
%evaluating a set of basis coordinates "tb" over a single "edge"

Pe = zeros(2, size(tb,2));

for ke = 1:size(edge,1)   
    Pe = Pe + bsxfun(@times, P(:, edge(ke,:)), C(ke,:)*tb);
end

end

function interpolants = internal_highOrderInt( PM, PS, bnd_edges_master, slavePoints, PsFun, Norder )
interpolants = [];

error('Most certainly does not work');

P_slave = PsFun( PS(:, slavePoints) );

%number of interpolation points
Nnext = ceil( (Norder+1) / 2 ); %no of points in positive direction
Nprev = floor( (Norder+1) / 2 ); %no of points in the negative direction

Nmaster = size(bnd_edges_master, 2);

%ordering edges
sorted_edges_master = order_edges(bnd_edges_master);
edgeLengths = sqrt( sum( (PM(:,sorted_edges_master(2,:)) - PM(:,sorted_edges_master(1,:))).^2, 1) );

%closest master edge to each slave point
closest_s2m = closestSegmentToNode( P_slave, ...
   PM(:, bnd_edges_master(1,:)),PM(:, bnd_edges_master(2,:)) );

Npoints = Nnext + Nprev;
Np_slave = size(closest_s2m, 2);

I = zeros(Np_slave*Npoints, 1); J = I; E = I;
ri = 1;
for kp_slave = 1:Np_slave
    me_init = closest_s2m(1, kp_slave); %initial master edge
       
    m_start = PM(:, bnd_edges_master(1,me_init)); %immediate prev master node
    m_end = PM(:, bnd_edges_master(2,me_init)); %immediate nexy master node

    l_point = norm( closest_s2m(2:3, kp_slave) - m_start );
    
    nodes_next = mod(((me_init+1):(me_init+Nnext-1)) - 1 , Nmaster) + 1;
    nodes_prev = mod(((me_init-0):-1:(me_init-Nprev)) - 1 , Nmaster) + 1;
    nodes_prev = fliplr(nodes_prev);
    
    nodes_all = [nodes_prev nodes_next];
    els = cumsum(edgeLengths(nodes_all)) - sum(edgeLengths(nodes_prev));
    
    %calculating interpolation weights (Lagrance interpolation polynomials
    ps = zeros(1, Nprev+Nnext);
    for kj = 1:(Nprev+Nnext)
        p = 1;
        for m = 1:(Nprev+Nnext)
            if m ~= kj
                p = p * (l_point - els(m)) / (els(kj) - els(m));
                if abs(els(kj) - els(m)) < 1e-4
                    warning('Suspiciously small angle. Check for dublicate nodes');
                end
            end
        end
        ps(kj) = p;
    end
    
    ps
    
    I(ri:(ri+Npoints-1)) = slavePoints(kp_slave) * ones(Npoints, 1);
    J(ri:(ri+Npoints-1)) = bnd_edges_master(1, nodes_all);
    E(ri:(ri+Npoints-1)) = ps;
    
    ri = ri + Npoints;
end

interpolants = [I J E];

end