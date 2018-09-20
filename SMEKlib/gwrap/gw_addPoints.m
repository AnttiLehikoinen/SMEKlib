function this = gw_addPoints(this, p, varargin)

Np = size(p,2);

%size check
if size(this.p, 2) < (this.N_points + Np)
    n_add = max(2*size(this.p,2), Np);
    this.p = [this.p zeros(3, n_add)];
end

%adding points
this.p(1:size(p,1), (this.N_points+1):(this.N_points+Np)) = p;
this.N_points = this.N_points + Np;
    
end