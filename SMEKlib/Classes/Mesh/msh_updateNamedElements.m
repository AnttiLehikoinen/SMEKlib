function msh = msh_updateNamedElements(msh, el_inds_new)
%msh_updateNamedElements
%
% (c) 2018 Antti Lehikoinen / Aalto University

keys = msh.namedElements.keys;
for k = 1:numel(keys)
    key = keys{k};
    val = msh.namedElements.get(key);
    if isa(val, 'cell')
        val2 = cell(size(val));
        for k2 = 1:numel(val)
            val2{k2} = el_inds_new( val{k2} );
        end
        val = val2;
    else
        val = el_inds_new(val);
    end
    msh.namedElements.set(key, val);
end

end