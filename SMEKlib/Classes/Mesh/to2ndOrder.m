function [p, t] = to2ndOrder(p, t, varargin)
%to2ndOrder transforms triangulation to second order.
% 
% [p, t] = to2ndOrder(p, t)
% [p, t] = to2ndOrder(p, t, edges, t2e)
% 
% (c) 2017 Antti Lehikoinen / Aalto University

if numel(varargin)
    edges = varargin{1};
    t2e = varargin{2};
else
    [edges, ~, t2e] = getEdges(t);
end

pnew = 0.5*p(:,edges(1,:)) + 0.5*p(:,edges(2,:));

Np = size(p,2);

p = [p pnew];
t = [t; Np+t2e];

end