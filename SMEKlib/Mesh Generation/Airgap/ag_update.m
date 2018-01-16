function [t_moving, p_virt] = ag_update(this, rotorAngle, varargin)

p_virt = this.p_virt;
rm = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)];
p_virt(:, this.n_moving) = rm*this.p_virt(:, this.n_moving);

%shifting indices in triangulation
%%{
nodeShift = -floor( (rotorAngle - this.shiftTol/2) / this.shiftTol ) - 1;
newPositions = mod(this.original_positions - 1 + nodeShift, numel(this.n_bnd) ) + 1;

t_moving = this.t_moving;
t_moving(this.inds_r) = this.n_bnd(newPositions);
return
%}

N_ag_r = 397;
edges = bsxfun(@minus, p_virt(:,this.n_bnd), p_virt(:,1) );
[~, nb_start] = min( sum(edges.^2,1) );
nb_new = this.n_bnd( mod((1:N_ag_r) + nb_start - 2, numel(this.n_bnd)) + 1);
%nb_new = this.n_bnd(1:N_ag_r);

N_static = size(this.p_virt, 2) - numel(this.n_moving);
msh_temp = struct('p', p_virt);
 
t_moving = singleLayerAGtriangulation_2(msh_temp, 1:N_static, nb_new);

end