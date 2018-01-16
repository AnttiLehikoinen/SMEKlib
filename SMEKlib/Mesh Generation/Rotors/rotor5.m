function msh = rotor5(msh,dim)

%Init the first sector
msh = calculate_rotor5(msh,dim);

%Boundaries
FL = [1 42 41 40 39 38 37 36];
LL = [1 24 25 26 27 28 29 30];
ag = [36 35 34 33 32 31 30];

%Copy the first sector
Nsec = dim.Qr/dim.num2;
[p, t, LLnew, agnew, ~] = replicate_sector_fixed(msh.p', msh.t', Nsec, dim.angleR(1), FL, LL, ag, 1);

msh.matel = repmat(msh.matel, Nsec, 1);

%finalization
msh.p = p';
msh.t = t';

msh.RC = reshape(find(msh.matel == 9999), [], Nsec)';

msh.index_p = size(msh.p,1);
msh.index_t = size(msh.t,1);
msh.n_ag_r = agnew';
msh.FL = setdiff(FL, 1);
msh.LL = LLnew;

end