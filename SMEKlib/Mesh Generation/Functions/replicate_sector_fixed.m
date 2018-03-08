function [p, t, n_master_final, n_ag, n_dir] = ...
    replicate_sector_fixed(p_sec, t_sec, N, theta, n_slave, n_master, n_ag, n_dir)
%replicate_sector replicates an elementary mesh sector.
%
% Call syntax
% [p, t] = replicate_sector(p_sec, t_sec, N, theta, n_slave, n_master)
% where
% p_sec and t_sec specify the elementary sector mesh, that is then replicated N
% times with theta radians rotation each. The nodes on the slave
% (clockwise) and master (counter-clockwise) boundaries of the elementary
% mesh are listed in n_slave and n_master.
%
% NOTE: If the segment contains an origin (0,0) node, it has to be listed
% in both n_slave and n_master;
% 
% (c) 2017 Antti Lehikoinen / Aalto University


%initializing nodes
Np_el = size(p_sec, 2); %number of nodes in the elementary mesh
n_rep = setdiff(1:Np_el, n_slave); %replicable nodes
Np_rep = numel(n_rep);
p = [p_sec zeros(2, (N-1)*Np_rep)];

%center node given?
n_center = intersect(n_slave, n_master);

if isempty(n_center)
    n_center = 0;
end

%initializing elements
Ne_el = size(t_sec, 2);
t = [t_sec zeros(3, (N-1)*Ne_el)];

%determining a replicable segment mesh
t_sec_rep = t_sec;
[inds_t_slave, locb] = ismember(t_sec, n_slave); %linear indices of n_slave within t_sec
inds_t_slave = find(inds_t_slave);
locb = locb(locb>0); %removing zero entries

ind_t_center = find( t_sec == n_center ); %linear index of center node, if any

%fixing indexing
pinds_new = zeros(1, Np_el);
pinds_new(n_rep) = 1:Np_rep;
t_sec_rep = reshape( pinds_new(t_sec_rep(:)), size(t_sec,1), []);

%replicating nodes and elements
for k = 2:N
   ra = (k-1)*theta; %rotation angle
   rM = [cos(ra) -sin(ra); sin(ra) cos(ra)]; %rotation matrix
   
   %replicating nodes
   p(:, (1:Np_rep) + Np_el + (k-2)*Np_rep) = rM*p_sec(:, n_rep);
   
   %replicating elements
   t_segment = t_sec_rep + (k-2)*Np_rep + Np_el;
   
   if k == 2
       t_segment(inds_t_slave) = n_master(locb);
   else
       t_segment(inds_t_slave) = pinds_new(n_master(locb)) + (k-3)*Np_rep + Np_el;
   end
   t_segment(ind_t_center) = n_center;
   t(:, (1:Ne_el) + (k-1)*Ne_el) = t_segment;
end

%updating nodal indices
n_master_final = n_master + Np_el + (N-2)*Np_rep;

n_ag_rep = pinds_new( n_ag(2:end) ) + Np_el;
n_ag = [n_ag repmat(n_ag_rep, 1, N-1)+kron(0:(N-2), Np_rep*ones(1,numel(n_ag_rep)))];

n_dir_rep = pinds_new( n_dir(2:end) ) + Np_el;
n_dir = [n_dir repmat(n_dir_rep, 1, N-1)+kron(0:(N-2), Np_rep*ones(1,numel(n_dir_rep)))];
    
end