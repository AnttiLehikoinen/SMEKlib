function this = gw_addLinePoints(this, xstart, xend, tol, varargin)

Np = ceil( norm(xend - xstart) / tol );
p = bsxfun(@plus, xstart(1:2,:), bsxfun(@times, xend(1:2,:)-xstart(1:2,:), linspace(0,1,Np)));

Np_orig = this.N_points;
this.addPoints(p);

Np = size(p, 2);
ldef = [1:(Np-1);2:Np] + Np_orig;
this.addLines(ldef, varargin{:});
end