classdef ExtrudedMachineMesh < ExtrudedPrismMesh
    properties
        
    end
    methods
        function msh3 = ExtrudedMachineMesh(msh2, zs)
            msh3 = msh3@ExtrudedPrismMesh(msh2, zs);
        end
    end
end