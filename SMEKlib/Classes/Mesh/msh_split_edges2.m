function msh = msh_split_edges2(msh, ed)
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
Nedges_new = numel(ed) + Ne_new;
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


%lazy for-loop implementation
ri_el = 1;
ri_edge_center = 1;
for kedge = 1:numel(ed)
    els = msh.e2t(:, ed(kedge)); %elements bordering this edge
    els = els( els > 0 );
    
    edge_old = ed( kedge ); %global edge index
    edge_new = Nedges + kedge; %index of the new edge
    
    n_new = Np + kedge; %index of new node
    
    %updating split edge
    n_1 = msh.edges(1, edge_old);
    n_2 = msh.edges(2, edge_old);
    
    msh.edges(:,edge_old) = sort([n_1; n_new]);
    msh.edges(:,edge_new) = sort([n_2; n_new]);
    
    for ke = 1:numel(els)
        %Splitting the element.
        % First, indices are determined in such a way that the old element
        % el is spanned by n_1, n_2, n_3, with n_1,n_2 belonging the old
        % edge to be split. Edges of the old element el are defined as
        % 1: n_1 --> n_2
        % 2: n_2 --> n_3
        % 3: n_3 --> n_1
        % Then, the element is split such that t(:,el) is replaced by
        % [n_2, n_3, n_new], n_new being the new edge-center node.
        % The new element [n_1; n_3; n_new] is appended to the list of
        % element.
        
        % determining lots of indices
        el = els(ke); %element index
        edgeind_element = find( msh.t2e(:, el) == edge_old ); %index of the edge to split, for this element
        el_new = Ne + ri_el; %index of new element
        edge_new_center = Nedges + numel(ed) + ri_edge_center; %index of the new edge splitting the element

        %determining ordering of the two other edges
        edges_in_element = mod( (edgeind_element + (0:2)) -1, 3) + 1;
        if ~intersect( n_2, msh.edges(:, msh.t2e(edges_in_element(2), el)) )
            edges_in_element([2 3]) = edges_in_element([3 2]);
        end        
        edges_in_element_global = msh.t2e(edges_in_element, el);

        %third corner node of the element
        n_3 = setdiff(msh.t(:, el), [n_1 n_2]);     

        %updating elements
        msh.t(:, el) = [n_2; n_3; n_new];
        msh.t(:, el_new) = [n_1; n_3; n_new];

        %updating element-splitting edge
        msh.edges(:,edge_new_center) = sort([n_3; n_new]);

        % updating element-to-edges
        %t2e_old = msh.t2e(:,el);
        t2e_old = edges_in_element_global;
        msh.t2e(:, el) = [t2e_old(2); edge_new_center; edge_new];
        msh.t2e(:, el_new) = [t2e_old(3); edge_new_center; edge_old];

        % updating edges-to-element:
        %element-splitting edge
        msh.e2t(:, edge_new_center) = sort([el; el_new]);

        %edge between n_1 and n_3
        temp = msh.e2t(:, edges_in_element_global(3)); %temp for clarity
        temp( temp==el ) = el_new; %updating the pointer to the old element
        msh.e2t(:, edges_in_element_global(3)) = sort(temp); %substituting back

        %edge between n_1 and n_new
        temp = msh.e2t(:, edge_old);
        msh.e2t(:, edge_new) = temp; %new edge shares the same elements as the old one did
        temp( temp==el ) = el_new;
        msh.e2t(:, edge_old) = sort(temp);
        

        %edge between n_1 and n_3 --> does not change
        %edge between n_1 and n_new --> does not change
        
        ri_el = ri_el + 1;
        ri_edge_center = ri_edge_center + 1;
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
    nn =  msh.namedElements.get(names{k});
    if iscell(nn)
        for kc = 1:numel(nn)
            temp( nn{kc} ) = 1;
            temp2 = temp( elementName_map );
            nn{kc} = find(temp2);
            temp = 0*temp;
        end
        msh.namedElements.set(names{k}, nn );
    else         
        temp( msh.namedElements.get(names{k}) ) = 1;
        temp2 = temp( elementName_map );
        msh.namedElements.set(names{k}, find(temp2) );
        temp = temp*0;
    end
end
msh.matel = msh.matel( elementName_map );

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