function S = closestSegmentToNode(P_nodes, P1, P2)
%given a set of node coordinates in P_nodes, and start and end coordinates
%of line segments P1-->P2, returns the closest segment S(1,:) to each node,
%and the coordinates S(2:3,:) of the closest point on segment

N_nodes = size(P_nodes, 2);
S = zeros(3, N_nodes);

tvectors = P2-P1; %tangent vectors of edges

%looping through points
for k_node = 1:N_nodes
    Cpoint = P_nodes(:, k_node); %point in question
    
    %First determine the closest point to Cpoint on ALL edges;
    %
    % express the unknown point as
    % p = P1 + kcoeff * tvector;
    % require that dot( Cpoint - p, tvector ) = 0;
    % solve kcoeff
    
    %closest point on lines determined by the edges
    kcoeffs = ( inner_dotProduct(Cpoint, tvectors) - dotProduct(P1, tvectors) ) ./ dotProduct(tvectors, tvectors);
    
    %restriction to the length of the edges
    kcoeffs = max(0, min(kcoeffs, 1));
    
    % Of these point candidates, determine the actually closest point
    Capp = P1 + bsxfun(@times, kcoeffs, tvectors);
    dists = sum( bsxfun(@plus, Cpoint, -Capp).^2, 1);
    
    [~,I] = min(dists);
    S(:, k_node) = [I; Capp(:,I)];
end

end

function E = inner_dotProduct(x, Y)
%dot product between a single vector x, and multiple vectors Y

E = x(1)*Y(1,:) + x(2)*Y(2,:);

end