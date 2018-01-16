function this = gw_addLineLoop(this, loops)

%size check
Nl = 1;
if size(this.ll, 2) < (this.N_lineloops + Nl)
    n_add = max(2*size(this.ll,2), Nl);
    this.ll = [this.l cell(1, n_add)];
end

%adding
this.ll{this.N_lineloops + 1} = loops;
this.N_lineloops = this.N_lineloops + Nl;

end