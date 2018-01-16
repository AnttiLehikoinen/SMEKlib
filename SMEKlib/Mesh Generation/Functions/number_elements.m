function [] = number_elements(p, t)
%number_elements plots element numbers.
% 
% number_elements(p, t) plots the indices of the elements (in t) in their
% centers (according to the nodal coordinates p.
%
% (c) 2017 Antti Lehikoinen

x = mean( reshape(p(1,t), size(t,1), []), 1);
y = mean( reshape(p(2,t), size(t,1), []), 1);

%plot(x, y, 'ko', 'Markersize', 12);
text(x, y, num2str((1:size(t,2))'), 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center','Fontsize',8);

end
