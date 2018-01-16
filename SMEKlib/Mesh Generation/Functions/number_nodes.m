function [] = number_nodes(p)
%number_nodes plots node numbers.
% 
% number_nodes(p) plots the indices of the nodes next to their coordinates
% in p.
%
% (c) 2017 Antti Lehikoinen / Aalto University

text(p(1,:), p(2,:), num2str((1:size(p,2))'), 'VerticalAlignment', 'bottom','Fontsize',8);

end
