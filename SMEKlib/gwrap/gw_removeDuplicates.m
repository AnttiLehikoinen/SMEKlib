function this = gw_removeDuplicates(this, varargin)
%gw_removeDuplicates remove duplicate points and lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removing duplicate nodes
if numel(varargin)
    TOL = varargin{1};
else
    TOL = 1e-4;
end

p = this.p(:, 1:this.N_points);

%p_id = round( p / TOL );
%[~, IA, IC] = unique(p_id', 'rows');

[~, IA, IC] = uniqueWithTol(p(1:2,:)', TOL, 'ByRows', true);

Np_new = numel(IA);
temp = 1:Np_new;

newInds = temp(IC);

this.p(:, 1:Np_new) = this.p(:, IA);
this.N_points = Np_new;

%propagating the new node indexing to the line definitions
this.l(:, 1:this.N_lines) = reshape( newInds( this.l(:,1:this.N_lines) ), 2, [] );

%FIXME can this change line directions?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removing duplicate lines
[this.l, linedirs] = sort(this.l, 1);
linedirs = linedirs(1,:);
linedirs( linedirs==2 ) = -1;

[~, IA, IC] = unique(this.l(:,1:this.N_lines)', 'rows');
Nl_new = numel(IA); temp = 1:Nl_new;
newInds = temp(IC);

this.l(:, 1:Nl_new) = this.l(:, IA);
this.N_lines = Nl_new;

%updating line loops
for kll = 1:this.N_lineloops
    ll = this.ll{kll};
    ll = newInds(abs(ll)) .* sign(ll) .* linedirs( abs(ll) );
    this.ll{kll} = ll;
end

%updating named lines if necessary
names = this.physicalLines.keys();
for k = 1:numel(names)
    ln = names{k};
    this.physicalLines.set(ln, newInds(this.physicalLines.get(ln)));
end

end
    