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
        misc
    end
    
    methods
        function this = AGtriangulation(msh, varargin)
            this = ag_constr(this, msh, varargin{:});
        end
        
        function Sag = get_AGmatrix(this, rotorAngle, varargin)
            Sag = ag_getAGmatrix(this, rotorAngle, varargin{:});
        end
        
        %CAUTION: changed output argument types
        %{
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
        %}
        function [tag_local, p, tag_global] = t_ag(this, varargin)
            disp('wfsd')
            if numel(varargin)
                rotorAngle = varargin{1};
            else
                rotorAngle = 0;
            end
            [t_mov, p] = this.update(rotorAngle);            
            tag_local = [this.t_const t_mov];
            if nargout == 0
                tag_global = reshape(this.agNodes_global(tag_local(:)), size(tag_local,1), []);
            else
                tag_global = [];
            end
        end
        
        function [t_moving, p_virt] = update(this, rotorAngle, varargin)
            [t_moving, p_virt] = ag_update(this, rotorAngle, varargin{:});
        end
        
        function setEccentricity(this, x)
            %sets eccentricity to the rotor
            % only works for two-band gap so far
            if ~isfield(this.misc, 'p_orig')
                %called for the first time --> initializing info
                this.misc.p_orig = this.p_virt;
                this.misc.agAngles_all = atan2( this.p_virt(2,:), this.p_virt(1,:) );
                this.misc.nr = setdiff(this.n_moving, this.n_bnd); %rotor boundary nodes
            end
            a = atan2(x(2), x(1)); %angle of the displacement
            nmid = this.n_bnd;
            %this.p_virt(:, nmid) = this.misc.p_orig(:,nmid) + ...
            %    r*[cos(this.misc.agAngles_all(nmid)-a); sin(this.misc.agAngles_all(nmid)-a)];
            this.p_virt(:, nmid) = this.misc.p_orig(:,nmid) + ...
                bsxfun(@times, x/2, abs(cos(this.misc.agAngles_all(nmid)-a)));
            nmid = this.misc.nr;
            this.p_virt(:, nmid) = bsxfun(@plus, this.misc.p_orig(:,nmid), x);
        end
    end
end
            