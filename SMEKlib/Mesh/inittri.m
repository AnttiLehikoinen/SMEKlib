function msh = inittri(p, t)
%inittri initializes a mesh struct
% 
% msh = inittri(p, t) initializes a mesh from nodes p and elements t
% The struct contains the following fields
%   p = nodes as 2xNp array
%   t = elements as 3xNe array
%   edges = edges as 2xNedges array
%   t2e = 3xNe array describing which edges each element has
%   e2t = 2xNedges array describing which elements each edge belongs to
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University,
%   based on the work of Antti Hannukainen and Mika Juntunen, 3.11.2005

if size(t,1)~=3
    t = t';
end

t = sort(t, 1);

%defining edges an' pals
[edges, e2t, t2e] = getEdges(t);

%struct
msh = struct('t', t, 'p', p, 'edges', edges, 'e2t', e2t, 't2e', t2e);

end