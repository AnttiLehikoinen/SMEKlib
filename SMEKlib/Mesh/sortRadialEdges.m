function n = sortRadialEdges(p, e)
%sortSegmentEdges sorts nodes of circumferentially.
%
%
% (c) 2018 Antti Lehikoinen / Smeklab

n = toRow(unique(e));

[~,I] = sort( sum(p(:,n).^2,1) );
n = n(I);

end