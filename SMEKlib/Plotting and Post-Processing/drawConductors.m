function [] = drawConductors(rs, Xc, varargin)
%drawConductors draws conductors as circles.
% 
% drawConductors(rs, Xc, args) draws round conductors with the radius rs
% each centered at columns of Xc. "args" contains the property-value pairs
% for the plotting function "rectangle".
% 
% Alternatively,
% drawConductors(rs, Xc, theta, args)
% rotates the conductors by theta radians.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if (numel(varargin) > 0) && isa(varargin{1}, 'double')
        rotA = varargin{1};
        Xc = [cos(rotA) -sin(rotA);sin(rotA) cos(rotA)] * Xc;
        plotArgs = {varargin{2:numel(varargin)}};
    else
        plotArgs = varargin;
end

for kc = 1:size(Xc,2)
    rectangle('Position', [Xc(:,kc)-rs; rs*[2;2]]', ...
        'Curvature', [1 1], plotArgs{:});
end

end