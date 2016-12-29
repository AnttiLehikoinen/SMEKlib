function S_struct = sparseEliminate(S_struct, columnwise, rowwise, el_table)
%sparseEliminate performs 1-on-1 elimination.
%
% Call syntax:
% S_struct = sparseEliminate(S_struct, columnwise, rowwise, el_table)
%
% A columnwise=true operation corresponds to setting
% S(:,el_table(1,:)) = S(:,el_table(2,:)) * diag(el_table(3,:))
%
% A rowwise=true operation corresponds to
% S(el_table(1,:),:) = diag(el_table(3,:)) * S(el_table(2,:))
%
% A faster implementation should be written
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

inds = 1:(S_struct.ri-1);

%{
S_struct.E(inds) = S_struct.E(inds) .* transpose(el_table(3,S_struct.I(inds)));
S_struct.I = el_table(2,S_struct.I(inds));

S_struct.E(inds) = S_struct.E(inds) .* transpose(el_table(3,S_struct.J(inds)));
S_struct.J = el_table(2,S_struct.J(inds));
return
%}
if rowwise
    [IA, IB] = ismember( S_struct.I(inds), el_table(1,:) ); 
    IA = IA(IB>0); IB = IB(IB>0);     
    S_struct.E(IA) = S_struct.E(IA) .* transpose(el_table(3,IB));
    S_struct.I(IA) = el_table(2,IB);
end
if columnwise
    [IA, IB] = ismember( S_struct.J(inds), el_table(1,:) );
    IA = IA(IB>0); IB = IB(IB>0);    
    S_struct.E(IA) = S_struct.E(IA) .* transpose(el_table(3,IB));
    S_struct.J(IA) = el_table(2,IB);
end

end    