function [] = mPlotG( G , varargin)
%mPlotG Plots geometry description matrix.
%   mPlotG(G, plotargs) plots the geometry description matrix G (used with 
%   Matlab's decsg function) with the plotargs arguments.
%   Functionality very limited for now.
%   
%   (c) 2017 Antti Lehikoinen / Aalto University

for k = 1:size(G,2)
    if (G(1,k)==2) || (G(1,k)==3)
        np = G(2,k);
        plot( G([3:(3+np-1) 3], k), G([(3+np):end 3+np], k), varargin{:});
    elseif G(1,k)==1
        r = G(4,k);
        rectangle('Position', [G(2,k)-r G(3,k)-r 2*r*[1 1]], 'Curvature', [1 1], ...
            'EdgeColor', varargin{1}, varargin{2:end});
    else
       warning('Geometry shape not yet implemented.')
    end

end

