function [] = edgePlot(P, EdgeSets, varargin)
%edgePlot plots edges.
% 
% edgePlot(P, EdgeSets, plotArgs) plots general-order edges. P contains the
% nodal coordinates, and EdgeSets contains the edge definitions on each
% column. EdgeSets can also be a cell array containing different-order
% edges.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if ~isa(EdgeSets, 'cell')
    EdgeSets = {EdgeSets};
end

Nplot = 40;
tplot = linspace(0, 1, Nplot);

Nsets = numel( EdgeSets );

for ke = 1:Nsets
    N_edgeOrder = size(EdgeSets{ke}, 1) - 1;
    
    C = internal_getRefShapeFunctions(N_edgeOrder);
    
    tb = bsxfun(@power, tplot, (0:N_edgeOrder)');
    
    [PX, PY] = internal_evaluateEdgeCoordinates(tb, C, P, EdgeSets{ke}, Nplot);
    
    plot(PX, PY, varargin{:});
end

end

function C = internal_getRefShapeFunctions(Norder)
%returns coefficients of the Lagrance shape functions of order Norder,
%evaluated at the reference edge

tpoints = linspace(0, 1, Norder + 1);

M = bsxfun(@power, tpoints', 0:Norder);


C = transpose( M \ eye(Norder+1) );

%valueMatrix = eye(Norder+1);
%C = transpose( M \ valueMatrix([1 Norder+1 2:Norder],:) );

end

function [PX, PY] = internal_evaluateEdgeCoordinates(tb, C, P, edges, Nplot)

PX = zeros(Nplot, size(edges,2));
PY = PX;

for ke = 1:size(edges,1)
    PX = PX + bsxfun(@times, P(1, edges(ke,:)), (C(ke,:)*tb)');
    PY = PY + bsxfun(@times, P(2, edges(ke,:)), (C(ke,:)*tb)');
end

end