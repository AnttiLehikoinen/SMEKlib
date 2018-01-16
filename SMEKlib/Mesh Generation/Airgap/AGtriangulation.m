classdef AGtriangulation < handle
    properties
        agNodes_global
        t_const
        t_moving
        p_virt
        msh_ag
        n_moving
        n_bnd
        S_const
        inds_r
        shiftTol
        original_positions
        el_table
    end
    
    methods
        function this = AGtriangulation(msh, varargin)
            this = ag_constr(this, msh, varargin{:});
        end
        
        function Sag = get_AGmatrix(this, rotorAngle, varargin)
            Sag = ag_getAGmatrix(this, rotorAngle, varargin{:});
        end
        
        function [tag, p, tag_local] = t_ag(this, varargin)
            if numel(varargin)
                rotorAngle = varargin{1};
            else
                rotorAngle = 0;
            end
            [t_mov, p] = this.update(rotorAngle);            
            tag_local = [this.t_const t_mov];
            if rotorAngle == 0
                tag = reshape(this.agNodes_global(tag_local(:)), size(tag_local,1), []);
            else
                tag = [];
            end
        end
        
        function [t_moving, p_virt] = update(this, rotorAngle, varargin)
            [t_moving, p_virt] = ag_update(this, rotorAngle, varargin{:});
        end
        
    end
end
            