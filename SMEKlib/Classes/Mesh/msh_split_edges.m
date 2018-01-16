function msh = msh_split_edges(msh, ed)
%split_edges Splits edges in two.
% 
% msh = split_edges(msh, ed)
% splits the edges specified in ed in two, and updates the element, node,
% and edge identities in msh.namedX
%
% NOTE: only one edge per element is split. In case two edges per element
% ke are given, only the edge with the lowest row-index in t2e(:,ke) is
% split.
% 
% The function also updates the listings in namedElements and namedEdges to
% correspond to the new mesh. The listings in namedNodes are updated in
% case a new node is added between two nodes in the listing.
%
% (c) 2017 Antti Lehikoinen / Aalto University

%determining new nodes
pnew = 0.5*msh.p(:, msh.edges(1,ed)) + 0.5*msh.p(:,msh.edges(2,ed));

%elements to split
t2split = [msh.e2t(1,ed) msh.e2t(2,ed)]; 
t2split = t2split( t2split>0 );

edges_split_orig = msh.edges(:, ed); %the split edges, expressed with nodes

%sizes
Ne_new = size(t2split, 2);
Nedges_new = Ne_new * 2;
Np_new = size(pnew, 2);
Np = size(msh.p,2); 
Ne = size(msh.t,2);
Nedges = size(msh.edges, 2);

%appending arrays
msh.p = [msh.p pnew];
msh.t = [msh.t msh.t(:, t2split)];
msh.edges = [msh.edges zeros(2, Nedges_new)];
msh.e2t = [msh.e2t zeros(2, Nedges_new)];
msh.t2e = [msh.t2e zeros(3, Ne_new)];

%relationship between the new edges and the new elements
[Lia, Locb] = ismember(msh.t2e(:, t2split), ed);

%lazy for-loop implementation
for ke = 1:size(t2split,2)
    % determining lots of indices
    el = t2split(:,ke); %element index
    
    edges2split_element = find(Lia(:, ke));
    for kee = 1:numel(edges2split_element)
        edgeind_element = edges2split_element(kee); %index of the edge to split, for this element

        el_new = ke + Ne;
        edgeind_list = Locb(edgeind_element, ke); %index of edge in the list of
            %edges to split

        edge_old = ed( edgeind_list ); %global edge index    
        edge_new = Nedges+edgeind_list; %index of the new edge
        edge_new_center = Nedges + Nedges_new/2 + edgeind_list; %index of the new edge splitting the element

        edges_in_element = mod( (edgeind_element + (0:2)) -1, 3) + 1; %ordering of edges
        edges_in_element_global = msh.t2e(edges_in_element, el);

        %nodes
        n_3 = setdiff(msh.t(:, el), msh.edges(:, edge_old));
        n_2 = intersect(msh.edges(:, edges_in_element_global(2)), msh.edges(:, edge_old));
        n_1 = intersect(msh.edges(:, edges_in_element_global(3)), msh.edges(:, edge_old));
        n_new = Np + edgeind_list;

        %updating elements
        msh.t(:, el) = [n_1; n_3; n_new];
        msh.t(:, el_new) = [n_2; n_3; n_new];

        %updating edges
        msh.edges(:,edge_old) = sort([n_1; n_new]);
        msh.edges(:,edge_new) = sort([n_2; n_new]);
        msh.edges(:,edge_new_center) = sort([n_3; n_new]);

        % updating element-to-edges
        t2e_old = msh.t2e(:,el);
        msh.t2e(:, el) = [t2e_old(2); edge_new_center; edge_old];
        msh.t2e(:, el_new) = [t2e_old(3); edge_new_center; edge_new];

        % updating edges-to-element:
        %element-splitting edge
        msh.e2t(:, edge_new_center) = sort([el; el_new]);

        %edge between n_2 and n_3
        temp = msh.e2t(:, edges_in_element_global(3)); %temp for clarity
        temp( temp==el ) = el_new; %updating the pointer to the old element
        msh.e2t(:, edges_in_element_global(3)) = sort(temp); %substituting back

        %edge between n_2 and n_new
        temp = msh.e2t(:, edge_old);
        temp( temp==el ) = el_new;
        msh.e2t(:, edge_new) = sort(temp);

        %edge between n_1 and n_3 --> does not change
        %edge between n_1 and n_new --> does not change
    end
end
msh.edges = sort(msh.edges, 1);
msh.e2t(:, ~msh.e2t(1,:)) = flipud( msh.e2t(:, ~msh.e2t(1,:)) );

% updating named elements and edges
elementName_map = [(1:Ne) t2split];
edgeName_map = [(1:Nedges) ed];

%elements:
names = msh.namedElements.keys();
temp = zeros(1, Ne);
for k = 1:numel(names)
    temp( msh.namedElements.get(names{k}) ) = 1;
    temp2 = temp( elementName_map );
    msh.namedElements.set(names{k}, find(temp2) );
    temp = temp*0;
end

%edges
names = msh.namedEdges.keys();
temp = zeros(1, Nedges);
for k = 1:numel(names)
    temp( msh.namedEdges.get(names{k}) ) = 1;
    temp2 = temp( edgeName_map );
    msh.namedEdges.set(names{k}, find(temp2) );
    temp = temp*0;
end

%nodes
%names = msh.namedNodes.keys();
names = {'Dirichlet'};
for k = 1:numel(names)
    nodes = msh.namedNodes.get(names{k});
    
    inds = find( ismember(edges_split_orig(1,:), nodes).*ismember(edges_split_orig(2,:), nodes) );
    msh.namedNodes.set(names{k}, ...
        [nodes Np+inds]);
end

end