classdef SurfaceMesh < handle
    properties
        p, t
    end
    
    methods
        function this = SurfaceMesh(p, t)
            this.p = p;
            this.t = t;
        end
        
        
        function [F, F0] = getMappingMatrix(this, varargin)
            if numel(varargin)
                els = varargin{1};
            else
                els = 1:size(this.t,2);
            end
            
            Fxy = [this.p(:, this.t(2,els))-this.p(:, this.t(1,els));
                this.p(:, this.t(3,els))-this.p(:, this.t(1,els))];
            F0 = this.p(:, this.t(1,els));
            
            Fz = crossProduct(Fxy(1:3, :), Fxy(4:6,:));
            Fz = bsxfun(@rdivide, Fz, sum(Fz.^2,1).^0.5);
            
            F = [Fxy; Fz];
        end
        
        function h = patch(this, els, varargin)
            if isempty(els)
                els = 1:size(this.t,2);
            end
            %reshape(this.p(1, this.t(:,els)), 3, [])
            %h = patch( reshape(this.p(1, this.t(:,els)), 3, []), ...
            %    reshape(this.p(2, this.t(:,els)), 3, []), ...
            %    reshape(this.p(3, this.t(:,els)), 3, []) );
            
            h = patch('Faces',this.t(:,els)', 'Vertices',this.p', varargin{:});
        end
        
        function [xq, wq] = quadpoints(~, integrand_order)
            [xq, wq] = get_2DtriangleIntegrationPoints(integrand_order);
            xq = [xq; zeros(1, numel(wq))];
        end
    end
    
end