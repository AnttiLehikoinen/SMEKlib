function n = sortSegmentEdges(p, e, varargin)
%sortSegmentEdges sorts nodes of circumferentially.
%
%
% (c) 2018 Antti Lehikoinen / Smeklab

n = toRow(unique(e));

[~,I] = sort( atan2(p(2,n), p(1,n)) );
n = n(I);

%sorting to begin from ~zero angle
if numel(varargin) && varargin{1}
    angles = atan2(p(2,n), p(1,n));
    ind_first = find( angles >= 0, 1 );
    n = n( mod( (ind_first:(ind_first+numel(n))-1) - 1, numel(n)) + 1);
    %n = n( mod( (ind_first:(ind_first+numel(n)))-1, numel(n)) + 1);
end
end