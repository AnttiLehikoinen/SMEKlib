function drawWindingConfiguration(W, dims)
%drawWindingConfiguration plots the winding configuration.
% 
% drawWindingConfiguration(W, dims) draws the winding configuration W for a
% machine with the dimensions dims (only Qs, q, and a needed).
% 
% NOTE: does not work very consistently for multi-pole machines.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

r_plot = 1;
r_mid = 2.5;

%virtual angle shift between layers
r_layerShift = 0.1;

angleShift = 0.5*pi/dims.Qs / size(W,1);


r_mid = r_mid - r_layerShift;
r_plot = r_plot - r_layerShift;

%a = dims.a;
a = max(abs(W(:))) / 3;
slotAngles = linspace(0/dims.Qs, (2-2/dims.Qs)*pi, dims.Qs);

phaseColors = {'k', 'r', 'b'};
colorFuns = {@(k_path)( 0.2*(k_path-1)*[1 1 1] );
    @(k_path)( [1 0 0] - (k_path-1)*[0.1 0 0] );
    @(k_path)( [0 0 1] - (k_path-1)*[0 0 0.1] )};

X = cell(3, 1); [X{:}] = deal(zeros(3, dims.q*a));
Y = cell(3, 1); [Y{:}] = deal(zeros(3, dims.q*a));

hold on;
for k_phase = 1:3
    ri = 1;
    for k_path = 1:a
        %[I_start,J_start] = internal_find( W ((k_phase-1)*a+k_path) );
        %[I_end, J_end] = internal_find( W == -((k_phase-1)*a+k_path) );
        
        k_ind = (-1)^(k_path) * ((k_phase-1)*a+k_path);
        [I_start, J_start, I_end, J_end] = internal_find(W, k_ind, dims);
        
        n = numel(I_start);
        
        
        angle_start = slotAngles(J_start) + (I_start-1)*angleShift;
        angle_end = slotAngles(J_end) + (I_end-1)*angleShift;        
        angle_mid = angle_start + 0.5*internal_angleDifference(angle_end, angle_start);
        
        X{k_phase}(1, ri:(ri+n-1)) = (r_plot + r_layerShift*I_start) .* cos(angle_start);
        Y{k_phase}(1, ri:(ri+n-1)) = (r_plot + r_layerShift*I_start) .* sin(angle_start);
        
        X{k_phase}(2, ri:(ri+n-1)) = (r_mid + r_layerShift*I_end) .* cos(angle_mid);
        Y{k_phase}(2, ri:(ri+n-1)) = (r_mid + r_layerShift*I_end) .* sin(angle_mid);
        
        X{k_phase}(3, ri:(ri+n-1)) = (r_plot + r_layerShift*I_end) .* cos(angle_end);
        Y{k_phase}(3, ri:(ri+n-1)) = (r_plot + r_layerShift*I_end) .* sin(angle_end);
        
        ri = ri + n;
    end
    plot( X{k_phase}, Y{k_phase}, 'Color', phaseColors{k_phase});
    %plot( X{k_phase}, Y{k_phase})
    plot( X{k_phase}([1 3],:), Y{k_phase}([1 3],:), 'Color', phaseColors{k_phase}, 'LineStyle', 'none', 'Marker', 'o');
end


end

function [I_start, J_start, I_end, J_end] = internal_find(W, k_path, dims)
%internal_find : custom find function;
% otherwise phases crossing from Q4 to Q1 would be plotted incorrectly

[~, inds] = find(W == k_path);

%quadrant-crossing detection
ind_start = 1;
for k = 2:size(W, 2)
    if any( ismember(W(:,k), k_path) ) && ~any( ismember(W(:,k-1), k_path) )
    %if (W(1,k) == k_path) && (W(1,k-1)~=k_path)
        ind_start = max(ind_start, k);
    end
end

if (max(inds)-min(inds)) > dims.q
    %ind_start = find( W(:,(dims.q+1):end) == k_path, 1) + dims.q;
else
    %ind_start = inds(1);
end

N = size(W,2);

%translating winding configuration to begin in the first slot of the path
W_temp = W(:, mod( (ind_start:(ind_start+N-1))-1, N)+1);

%finding slots and unrolling
[I_start, J_start] = find(W_temp == k_path);
I_start = internal_2row( I_start );
J_start = internal_2row( mod(J_start + ind_start, N) + 1 );

[I_end, J_end] = find(W_temp == -k_path);
I_end = internal_2row( I_end );
J_end = internal_2row( mod(J_end + ind_start, N) + 1 );

end

function x = internal_2row(x)
if iscolumn(x)
    x = x';
end
end



function d = internal_angleDifference(a1, a2)

%calculates the signed difference a1-a2 between the angles
d = mod(a1 - a2 + pi,2*pi*(1+1e-9)) - pi;

end