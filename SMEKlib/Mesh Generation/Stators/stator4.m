function msh = stator4(msh,dim)
A = size(find(msh.matel == dim.SSM1));

%Init the first sector
[p,t,m] = calculate_stator4(dim,msh);

%Boundaries
FL = [32 31 30 29 28];
LL = [22 23 24 25 26];
cir = [28 27 26];
ag = [32 18 19 2 1 12 20 21 22];

Np_old = size(msh.p, 1);

Nsec = dim.Qs/dim.num2;
[p, t, LLnew, agnew, cirnew] = replicate_sector_fixed(p', t', Nsec, dim.angleS(1), FL, LL, ag, cir);

%Finalize
msh.matel = [msh.matel;repmat(m, Nsec, 1)];

if dim.SSM1 ~= dim.SSM2
    SC = find(msh.matel == dim.SSM1);
    SC = SC(A1+1:end);
    SC1 = find(msh.matel == dim.SSM2);
    SC1= SC1(A2+1:end);
    SC = [SC;SC1];
else
    SC = find(msh.matel == dim.SSM1);
    SC = SC(A+1:end);
end

%Finalize
msh.cir = [1 Np_old+cirnew];
msh.p = [msh.p;p'];
msh.t = [msh.t; Np_old + t'];
msh.n_ag_s = Np_old + agnew';
msh.FL = [msh.FL'; Np_old+FL'];
msh.index_p = size(msh.p,1);
msh.LL = [msh.LL'; LLnew'];
SC = sort(SC);
msh.SC = reshape(SC, [], Nsec)';

msh.matel(msh.matel == 9999) = dim.SWM;
end