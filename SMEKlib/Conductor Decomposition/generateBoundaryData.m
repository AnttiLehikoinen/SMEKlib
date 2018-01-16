function boundaryData = generateBoundaryData(PM, PS, bnd_edges_master, bnd_edges_slave, PsFun)

%ordering edges
sorted_edges_master = order_edges(bnd_edges_master);
sorted_edges_slave = order_edges(bnd_edges_slave);

%establishing relationship between sorted and original edges
%[~,L_master] = ismember(bnd_edges_master', sort(sorted_edges_master,1)', 'rows');
%[~,L_slave] = ismember(bnd_edges_slave', sort(sorted_edges_slave,1)', 'rows');
%disp(bnd_edges_master)
%disp(sort(sorted_edges_master(:,L),1))

[~, L_master] = ismember(sort(sorted_edges_master,1)', sort(bnd_edges_master,1)', 'rows');
[~, L_slave] = ismember(sort(sorted_edges_slave,1)', sort(bnd_edges_slave, 1)', 'rows');

%obtaining necessary node coordinates
P_master = PM(:, sorted_edges_master(1,:));
P_slave = PsFun( PS(:, sorted_edges_slave(1,:)) );

%checking orientation
%%{
if ispolycw(P_master(1,:), P_master(2,:))
    disp('Master edges not organized counter-clockwise. Flipping...')
    L_master = fliplr(L_master');
    sorted_edges_master = fliplr(sorted_edges_master([2 1],:));
    P_master = PM(:, sorted_edges_master(1,:));
end
if ispolycw(P_slave(1,:), P_slave(2,:))
    disp('Slave edges not organized counter-clockwise. Flipping...')
    L_slave = fliplr(L_slave');
    sorted_edges_slave = fliplr(sorted_edges_slave([2 1],:));
    P_slave = PsFun( PS(:, sorted_edges_slave(1,:)) );
end
%}


%changing references
[C_master, ~, IC_master] = unique( sorted_edges_master(:), 'stable' );
temp = 1:numel(C_master); 
sorted_edges_master = reshape(temp(IC_master), 2, []);

%sorted_edges_slave
[C_slave, ~, IC_slave] = unique( sorted_edges_slave(:), 'stable' );
temp = 1:numel(C_slave); 
sorted_edges_slave = reshape(temp(IC_slave), 2, []);
%reshape( C_slave(sorted_edges_slave(:)), 2, [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding closest segment to each point

%closest slave edge to each master point
closest_m2s = closestSegmentToNode( P_master(:, sorted_edges_master(1,:)), ...
    P_slave(:, sorted_edges_slave(1,:)), P_slave(:, sorted_edges_slave(2,:)) );

%closest master edge to each slave point
closest_s2m = closestSegmentToNode( P_slave(:, sorted_edges_slave(1,:)), ...
   P_master(:, sorted_edges_master(1,:)),P_master(:, sorted_edges_master(2,:)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining the overlapping segments node-by-node
tol = 1e-8; %tolerance for node distance

Nmaster = size(P_master, 2);
Nslave = size(P_slave,2);

Nmax = 2*max( Nmaster, Nslave );

temp_data = zeros(2, Nmax);
segment_points = zeros(size(temp_data));
ri = 1;
for k = 1:(Nslave+1)
    slaveNode_start = mod(k-1,Nslave) + 1;

    slaveNode_end = mod(k+1-1,Nslave) + 1;

    temp_data(:,ri) = [slaveNode_start; 2]; 
    segment_points(:,ri) = P_slave(:, slaveNode_start);
    ri = ri + 1;

    %master nodes falling on this interval
    inds = find( closest_m2s(1,:) == slaveNode_start);
    dists_start = sum( bsxfun(@plus, closest_m2s(2:3,inds), -P_slave(:,slaveNode_start)).^2 ).^0.5;
    dists_end = sum( bsxfun(@plus, closest_m2s(2:3,inds), -P_slave(:,slaveNode_end)).^2 ).^0.5;
    inds = inds( (dists_start>tol) & (dists_end>tol) );    
    if any(inds)
        Nincr = numel(inds);    
        temp_data(:,ri:(ri+Nincr-1)) = [inds; ones(1,Nincr)]; 
        segment_points(:,ri:(ri+Nincr-1)) = P_master(:,inds);
        ri = ri + Nincr;
    end
end
temp_data = temp_data(:, 1:(ri-1) );
segment_points = segment_points(:, 1:(ri-1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating boundary data
Nsegs = size(temp_data, 2)-1;

P_centers = ( segment_points(:,2:end)+segment_points(:,1:(end-1)) )/2; %centers of segments
%closest slave edge to each segment center
closest_SE = closestSegmentToNode(P_centers , ...
    P_slave(:, sorted_edges_slave(1,:)), P_slave(:, sorted_edges_slave(2,:)) );

%closest master edge to each segment center
closest_ME = closestSegmentToNode(P_centers , ...
    P_master(:, sorted_edges_master(1,:)), P_master(:, sorted_edges_master(2,:)) );

%finally generating the boundary data
boundaryData = zeros(10, Nsegs);
for k = 1:Nsegs
    startNode = temp_data(1, k);
    nextColumn = mod(k+1-1,ri) + 1;
    nextNode = temp_data(1, nextColumn);

    %dealing with start node
    if temp_data(2, k) == 2
        %boundaryData(1:3, k) = closest_s2m(:, startNode);
        boundaryData(1:3, k) = [closest_ME(1,k); closest_s2m(2:3, startNode)];
        %boundaryData(1:3, k) = [closest_ME(1,startNode); P_slave(:, startNode)];
        boundaryData(6:8, k) = [startNode; P_slave(:, startNode)];
    else
        boundaryData(1:3, k) = [startNode; P_master(:,startNode)];
        %boundaryData(6:8, k) = closest_m2s(:, startNode);
        boundaryData(6:8, k) = [closest_SE(1,k); closest_m2s(2:3, startNode)];
        %boundaryData(6:8, k) = [closest_SE(1,startNode); P_master(:,startNode)];
    end

    %dealing with end node
    if temp_data(2, nextColumn) == 2
        boundaryData(4:5, k) = closest_s2m(2:3, nextNode);
        boundaryData(9:10, k) = P_slave(:,nextNode);
    else
        boundaryData(4:5, k) = P_master(:,nextNode);
        boundaryData(9:10, k) = closest_m2s(2:3, nextNode);
    end
end

boundaryData(1, :) = L_master( boundaryData(1, :) );
boundaryData(6, :) = L_slave( boundaryData(6, :) );

end
