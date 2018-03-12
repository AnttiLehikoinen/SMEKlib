function n = sortSegmentEdges(p, e)
%sortSegmentEdges sorts nodes of segment
%
%
% (c) 2018 Antti Lehikoinen / Smeklab

n = toRow(unique(e));

[~,I] = sort( atan2(p(2,n), p(1,n)) );
n = n(I);

end