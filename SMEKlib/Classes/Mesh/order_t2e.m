function t2e_ordered = order_t2e(msh)
%order_t2e orders the triangles-to-edges matrix.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

t2e_ordered = msh.t2e;

for kn = 1:3
    inds = find( msh.t(kn,:) ~= msh.edges(1, msh.t2e(kn,:)) );
    t2e_ordered(kn, inds) = -t2e_ordered(kn, inds);
end

end