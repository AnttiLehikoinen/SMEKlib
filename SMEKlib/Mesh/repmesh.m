function [pnew, tnew, varargout] = repmesh(p, t, Pfun, varargin)
%repmesh nightmare function for extending meshes. Work-in-progress.
% 
% [pnew, tnew, otherStuff] = repmesh(p, t, Pfun, stuff) tries to
% replicate the initial mesh based on the node-moving function Pfun an'
% stuff, and hopefully parse the triangulations together on the boundary.
%
% Call syntax 1:
% [pnew, tnew] = repmesh(p, t, Pfun, inds_p_master, inds_p_slave, N_rep)
%   Pfun = a function handle for moving the nodes (e.g. translation)
%   inds_p_master = which nodes of the original triangulation are preserved on
%      the boundary where the original and replicated triangulations meet
%   inds_p_slave = which nodes are slaved to the aforementioned master
%      nodes
%
% Additionally, call with
% repmesh(p, t, Pfun, inds_p_master, inds_p_slave, N_rep, 't', t1, 't',
% t2, 'p', p1, 'p', p2) to get the indices of some elements t1,t2 and nodes
% p1,p2 in the new triangulation. In this case, the output is
% [pnew, tnew, t1_new, t2_new, p1_new, p2_new]
%
% Call syntax 2:
% repmesh(p, t, N_sectors) to replicate a sector mesh for the entire
% cross-section. Master and slave nodes are determined automatically,
% fortune permitting. Element and node indices can also be supplied if
% necessary.
%
% Copyright (c) 2015-2016 Antti Lehikoinen / Aalto University


if isa(Pfun, 'function_handle')
    wrapMesh = false;
    inds_p_master = varargin{1};
    inds_p_slave = varargin{2};
    N_rep = varargin{3};
    inds_list = varargin{4:end};
else
    wrapMesh = true;
    N_rep = Pfun - 1;
    angle_sector = 2*pi/(N_rep + 1);
    Pfun = @(X)( [cos(angle_sector) -sin(angle_sector);sin(angle_sector) cos(angle_sector)]*X);
    inds_list = varargin;
    
    p_angles = atan2(p(2,:), p(1,:));
    p_angles = p_angles - min(p_angles); %rotating sector to the first quadrant
    
    tol = 1e-3;
    ind_p_center = find( sum(p.^2,1).^0.5 < 0.5e-3 );
    inds_p_slave = find( abs(p_angles) < tol );
    inds_p_slave = setdiff(inds_p_slave, ind_p_center); %possible center node does not change
    
    inds_p_master = find( abs(p_angles - angle_sector) < tol );
end


%figure(1); clf; hold on;
%plot(p(1, inds_p_slave), p(2,inds_p_slave), 'ko');
%plot(p(1, inds_p_master), p(2,inds_p_master), 'rx');

%repeats mesh
Np_old = size(p, 2);
Ne_old = size(t, 2);

pnew = zeros(2, (N_rep+1)*Np_old);
tnew = zeros(3, (N_rep+1)*Ne_old);

ptemp = p; pnew(:, 1:Np_old) = ptemp;
tnew(:, 1:Ne_old) = t;
for krep = 1:(N_rep)
    ptemp = Pfun(ptemp);
    pnew(:, (krep*Np_old+1):((krep+1)*Np_old)) = ptemp;
    tnew(:, (krep*Ne_old+1):((krep+1)*Ne_old)) = krep*Np_old + t;
end

Np_new = size(pnew,2);

temp_pinds = 1:Np_new;

inds_p2kill = kron(1:N_rep, Np_old * ones(1, numel(inds_p_slave))) + repmat(inds_p_slave, 1, N_rep);
temp_pinds(inds_p2kill) = kron(0:(N_rep-1), Np_old * ones(1, numel(inds_p_master))) + ...
    repmat(inds_p_master, 1, N_rep);

%dealing with the possible center node
if any(ind_p_center)
    temp_pinds(ind_p_center + (1:N_rep)*Np_old) = ind_p_center;
end

if wrapMesh
    temp_pinds(inds_p_slave) = (N_rep*Np_old) + inds_p_master;
end
    

[p2keep, IA, IC] = unique(temp_pinds);

temp = 1:numel(IA); newInds = temp(IC);

pnew = pnew(:, p2keep);
tnew = reshape( newInds(tnew(:)), 3, []);


ri = 1; index_type = 1;
for kout = 1:numel(inds_list)
    if isa(inds_list{kout}, 'char') && strcmp(inds_list{kout}, 't')
        index_type = 1;
        continue;
    elseif isa(inds_list{kout}, 'char') && strcmp(inds_list{kout}, 'p')
        index_type = 2;
        continue;
    end
    
    if index_type == 1
        varargout{ri} = kron(0:N_rep, Ne_old*ones(1, numel(inds_list{kout}))) + ...
            repmat(inds_list{kout}, 1, N_rep+1);
        ri = ri + 1;
    else
       varargout{ri} = kron(0:N_rep, Np_old*ones(1, numel(inds_list{kout}))) + ...
            repmat(inds_list{kout}, 1, N_rep+1);
        ri = ri + 1;
    end
end
        


end