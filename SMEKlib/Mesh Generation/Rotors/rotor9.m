function msh = rotor9(msh,dim)
%This is still unfinished

%Init the first sector
[p,t,m,FL,LL,ag] = calculate_rotor9(msh,dim);
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
%{
%Init the first sector
msh = calculate_rotor9(msh,dim);

%Boundaries
FL = [1 68 67 66 65 64 63 62];
LL = [1 46 47 48 49 50 51 52];
ag = [62 61 60 59 58 57 56 55 54 53 52]';
t_len = numel(msh.t)/3+1;
ag_len = numel(ag);

%Copy the first sector
[p, t] = replicate_sector(msh.p', msh.t', dim.Qr/dim.num, dim.angleR(1), FL, LL);
msh.p = p';
msh.t = t';

%Make matel & air-gap vectors
[m,ag] = ultimate_vector_expander(msh.matel,ag(2:end),dim.Qr/dim.num,msh.index_p-numel(LL),[]);
msh.matel = m;
msh.RC = vec2mat(find(msh.matel == 9999),30);
msh.matel(msh.matel == 9999) = dim.RSM;

%Index corrections
[msh.t,ag] = awesome_index_corrector(0,msh,[],ag,[],t_len,ag_len,LL,0);

%Finalize mesh
msh.index_p = size(msh.p,1);
msh.index_t = size(msh.t,1);
msh.n_ag_r = [62 ;ag];
msh.FL = [68 67 66 65 64 63 62];
msh.LL = [msh.index_p-15 msh.index_p-14 msh.index_p-13 msh.index_p-12 msh.index_p-11 msh.index_p-10 msh.index_p-9];
%}
end