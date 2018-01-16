classdef MeshBase3D < handle
    %MeshBase3D a base class for a 3D mesh.
    % 
    % Subject to changes.
    %
    % (c) 2017 Antti Lehikoinen / Aalto University
    
    properties
        elementType
        nodes, elements
        edges, faces,
        faces2elements, elements2faces, elements2edges
        namedNodes, namedElements, namedFaces, info
    end
    
    methods
        function msh = MeshBase3D()
            msh.namedNodes = SLContainer();
            msh.namedElements = SLContainer();
            msh.namedFaces = SLContainer();
            msh.info = SLContainer();
        end
        
    end
    
end