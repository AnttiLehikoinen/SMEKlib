classdef SurfaceLineMesh < handle
    properties
        p, t
    end
    
    methods
        function this = SurfaceLineMesh(p,t)
            this.p = p;
            this.t = t;
        end
        
        function [F, F0] = getMappingMatrix(this, varargin)
            if numel(varargin)
                els = varargin{1};
            else
                els = 1:size(this.t,2);
            end
            F0 = this.p(:, this.t(1,els));
            
            F1 = this.p(:, this.t(2,els)) - F0;
            
            %getting normal vector
            F2 = zeros(3, numel(els));
            inds = els;
            %looping until all new normals are long enough
            while numel(inds) > 0
                %new random vector
                F2(:, inds) = crossProduct(F1(:,inds), rand(3, numel(inds)));
                
                inds = els( sum(F2(:,inds).^2,1) < 1e-10 );
            end
            %normalizing
            F2 = bsxfun(@rdivide, F2, sum(F2.^2,1).^0.5);
            
            %final vector
            F3 = crossProduct(F1, F2);
            F3 = bsxfun(@rdivide, F3, sum(F3.^2,1).^0.5);
            
            F = [F1; F2; F3];
        end
    
        function h = patch(this, els, varargin)
            if isempty(els)
                els = 1:size(this.t,2);
            end
            
            h = plot3( [this.p(1, this.t(1,els)); this.p(1, this.t(2,els))], ...
                [this.p(2, this.t(1,els)); this.p(2, this.t(2,els))], ...
                [this.p(3, this.t(1,els)); this.p(3, this.t(2,els))], ...
                varargin{:});
        end
        
        function [xq, wq] = quadpoints(~, integrand_order)
            [xq, wq] = get_1DQuadPoints_10(integrand_order);
            xq = [xq; zeros(2, numel(wq))];
        end
    end
    
    
end