function this = gw_addLines(this, lines)

%size check
Nl = size(lines, 2);
if size(this.l, 2) < (this.N_lines + Nl)
    n_add = max(2*size(this.l,2), Nl);
    this.l = [this.l zeros(2, n_add)];
end

%adding
this.l(:, (this.N_lines+1):(this.N_lines+Nl)) = lines;
this.N_lines = this.N_lines + Nl;
end
