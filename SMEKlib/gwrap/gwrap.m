classdef gwrap < handle
    properties
        gpath
        N_points, N_lines, N_lineloops, N_surfaces,
        p, l, ll
        surfaces, holesInSurfaces, physicalLines
    end
    
    methods
        function this = gwrap(gmsh_path)
            
            if gmsh_path(end)~='\'
                 this.gpath = [gmsh_path '\'];
            else
                this.gpath = gmsh_path;
            end
           
            
            %initializing numbers
            this.N_points = 0;
            this.N_lines = 0;
            this.N_lineloops = 0;
            this.N_surfaces = 0;
            
            %initializing arrays
            this.p = zeros(3, 1000);
            this.l = zeros(2, 1000);
            this.ll = cell(1, 1000);
            this.surfaces = SLContainer();
            this.holesInSurfaces = SLContainer();
            this.physicalLines = SLContainer();
        end
        
        function this = addPoints(this, p)
            this = gw_addPoints(this, p);
        end
        function this = addLines(this, l)
            this = gw_addLines(this, l);
        end
        function this = addLineloop(this, ll)
            this = gw_addLineLoop(this, ll);
        end
        function this = addSurface(this, points, surfacename, varargin)
            this = gw_addSurface(this, points, surfacename, varargin{:});
        end
        function this = addPcws(this, varargin)
            %Add surface piecewise. 
            %   help gw_addPcws for more details.
            this = gw_addPcws(this, varargin{:});
        end
        function this = addpcw(this, varargin)
            this = gw_addpcw(this, varargin{:});
        end
        
        function this = removeDuplicates(this, varargin)
            this = gw_removeDuplicates(this, varargin{:});
        end
        
        function this = writeFile(this, varargin)
            gw_writeFile(this, varargin{:});
        end
        function this = mesh(this, varargin)
            gw_mesh(this, varargin{:});
        end
        function [p, t, Surfaces] = loadMesh(this, varargin)
            [p, t, Surfaces] = gw_loadMesh(this, varargin{:});
        end
        
            
        
        function this = plotSurface(this, surfacename, varargin)
            gw_plotSurface(this, surfacename, varargin{:});
        end
        function this = fillSurface(this, surfacename, varargin)
            gw_fillSurface(this, surfacename, varargin{:});
        end
    end
end