function msh = mesh_stator(msh, dim, varargin)
if dim.SSM1 ~= dim.SSM2
    A1 = size(find(msh.matel == dim.SSM1));
    A2 = size(find(msh.matel == dim.SSM2));
else
    A = size(find(msh.matel == dim.SSM1));
end

%Init the first sector
if numel(varargin) && isa(varargin{1}, 'function_handle')
    [p,t,m,FL,LL,cir,ag] = varargin{1}(dim);
else
    [p,t,m,FL,LL,cir,ag] = calculate_stator1(dim);
end


Np_old = size(msh.p, 1);
Nsec = dim.Qs/dim.num2;

[p, t, LLnew, agnew, cirnew] = replicate_sector_fixed(p', t', Nsec, dim.angleS(1), FL, LL, ag, cir);

msh.matel = [msh.matel;repmat(m, Nsec, 1)];

if dim.SSM1 ~= dim.SSM2
    SC1 = find(msh.matel == dim.SSM1);
    SC2 = find(msh.matel == dim.SSM2);
    %SC = SC(A1+1:end);
    %SC1 = find(msh.matel == dim.SSM2);
    %SC1= SC1(A2+1:end);
    %SC = [SC; SC1];
    
    %SC = sort(SC);
    SC1 = reshape(SC1, [], Nsec);
    SC2 = reshape(SC2, [], Nsec);
    SC = cell(1, 2*Nsec);
    for k = 1:2*Nsec
        if mod(k,2)
            SC{k} = SC1(:, floor((k-1)/2)+1 )';
        else
            SC{k} = SC2(:, floor((k-1)/2)+1 )';
        end
    end
    
    msh.matel( msh.matel == dim.SSM2 ) = dim.SSM1;
else
    SC = find(msh.matel == dim.SSM1);
    SC = SC(A+1:end);
    
    %SC = sort(SC);
    SC = reshape(SC, [], Nsec);
end

%Finalize
%msh.cir = [1 Np_old+cirnew];
msh.cir = [msh.cir Np_old+cirnew];
msh.p = [msh.p;p'];
msh.t = [msh.t; Np_old + t'];
msh.n_ag_s = Np_old + agnew';
msh.FL = [msh.FL'; Np_old+FL'];
msh.index_p = size(msh.p,1);
msh.LL = [msh.LL'; LLnew'];
msh.SC = SC;
msh.matel(msh.matel == 9999) = dim.SWM;

msh.FL = setdiff(msh.FL, msh.cir);
msh.LL = setdiff(msh.LL, msh.cir);

end