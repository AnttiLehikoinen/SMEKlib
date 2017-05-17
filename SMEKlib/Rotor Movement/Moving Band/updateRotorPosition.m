function [t_ag, nodeShift, pnew] = updateRotorPosition(rotorAngle, msh)
%updateRotorPosition updates the air-gap mesh
% 
% [t_ag, nodeShift, pnew] = updateRotorPosition(rotorAngle, msh)
%
% NOTE: if msh has a field symmetrySectors, pnew contains air-gap nodes
% (real and virtual both) only, and t_ag is indexed to refer to them.
%
% Copyright (c) 2015-2016 Antti Lehikoinen / Aalto University

%number of nodes to "skip"
nodeShift = -floor( (rotorAngle - 1*msh.bandData.shiftTol/2) / msh.bandData.shiftTol ) - 1;

%shifting indices in the sorted list
newPositions = mod(msh.bandData.originalPositions_rotor - 1 + nodeShift, msh.bandData.N_ag_r ) + 1;

t_ag = msh.bandData.t_ag;
t_ag(msh.bandData.inds_r) = msh.bandData.sortedNodes_rotor(newPositions);

if ~isfield(msh, 'symmetrySectors')
    t_ag = reshape(msh.bandData.agNodes_global(t_ag(:)), 3, []); %shifting to global indexing
    
    pnew = msh.p;
    inds_p_rotor = msh.bandData.agNodes_global( msh.bandData.sortedNodes_rotor );
    pnew(:, inds_p_rotor) = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)] * ...
        pnew(:, inds_p_rotor);
else
    pnew = msh.bandData.p_ag_virt;
    pnew(:, msh.bandData.sortedNodes_rotor) = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)] * ...
        pnew(:, msh.bandData.sortedNodes_rotor);
end    

end