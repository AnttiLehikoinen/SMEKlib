function msh = rotor20(msh,dim)

%Init the first sector
[p,t,m,FL,LL,ag] = calculate_rotor20(dim);
Nsec = dim.Qr/dim.num;
[p, t, LLnew, agnew, ~] = replicate_sector_fixed(p', t', Nsec, dim.angleR(1), FL, LL, ag, []);
msh.matel = repmat(m, Nsec, 1);

%Finalize mesh
msh.p = p';
msh.t = t';
msh.index_p = size(msh.p,1);
msh.index_t = size(msh.t,1);
msh.n_ag_r = agnew';
msh.FL = FL(2:end);
msh.index_p = size(msh.p,1);
msh.LL = LLnew(2:end)-1;
msh.RC = reshape(find(msh.matel == 9999), [], Nsec);
msh.matel(msh.matel == 9999) = dim.RSM;
msh.matel(msh.matel == 999) = dim.RO;

end