function el = get_elementsInDomain(bt, minRegion)
%get_elementsInDomain gets elements in domain.
% 
% el = get_elementsInDomain(bt, elementRegion) returns the indices of
% elements in each domain. bt is the boolean table between minimal domains
% and domains (see Matlab's decsg), and minRegion is the minimal region
% each element belongs to.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

el = cell(1,size(bt,2));

for k = 1:size(bt,2)
    el{k} = find( ismember(minRegion, find(bt(:,k))) );
end
    

end