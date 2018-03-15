function [hp, ht] = pointplot(x, name, varargin)
%pointplot demonstrative plot of gwrap control point.
%
% [hp, ht] = pointplot(x, name, plot_args)
% 
% (c) 2018 Antti Lehikoinen / Smeklab

hp = plot(x(1,:), x(2,:), varargin{:});

%plotting point name
ht = text(x(1), x(2), name, 'VerticalAlignment', 'bottom','Fontsize',8);

end